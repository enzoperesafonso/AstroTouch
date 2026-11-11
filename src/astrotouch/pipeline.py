"""
The model-building pipeline, shared by the command line and the GUI:

    load -> remove stars / painted areas -> clean noise -> downsample + stretch
         -> compact islands -> flatten painted areas -> smooth -> markers
         -> texture -> mesh (border, label strip, label)

Hand-painted layers (paint-out, protect, texture) live on a "mask grid": the
input image at most MAX_MASK_GRID pixels across (raw pixels are grouped in
square cells of mask_scale() pixels), so painting stays fast on huge mosaics.
Cells are far finer than anything a finger can feel.
"""
from dataclasses import dataclass, field, asdict
import logging
import threading

import numpy as np
from scipy import ndimage

try:
    from . import core, shaping, sources
    from .clean import clean_image, beam_from_header
except ImportError:  # imported as a top-level module
    import core
    import shaping
    import sources
    from clean import clean_image, beam_from_header

logger = logging.getLogger(__name__)

MAX_MASK_GRID = 2048
STRETCHES = ("linear", "log", "asinh", "tactile")
PAINT_MODES = ("fill", "flatten")
LABEL_PADS = ("base", "raised")


@dataclass
class Settings:
    """Every option of a model. Defaults are shared by the CLI and the GUI."""
    # Size
    longest_side_mm: float = 150.0
    max_height_mm: float = 10.0
    base_thickness_mm: float = 2.0
    # Image
    clip: float = None                 # None: 1%, or 0 with cleaning or the tactile stretch
    smooth: float = 1.0                # Gaussian sigma, model pixels
    downsample: int = 1
    downsample_method: str = "mean"
    ridge: float = 0.0
    min_pixel_mm: float = 0.3
    nan_value: float = None
    stretch: str = "linear"
    soften: float = 0.0
    invert: bool = False
    # Tactile stretch
    edge: float = None
    edge_step: float = 0.2
    contrast: float = 1.5
    peak_percentile: float = 99.9
    # Noise cleaning
    clean: bool = False
    beam: tuple = None
    detect_sigma: float = 3.0
    peak_sigma: float = 5.0
    min_island: float = None
    smooth_scales: tuple = None
    noise_box: int = None
    # Stars and compact sources
    mask_stars: bool = False
    gaia_cache: str = None
    remove_points: bool = False
    point_box: int = 9
    point_sigma: float = 5.0
    max_point_area: float = 60.0
    masks: tuple = ()                  # (x, y, r) discs, raw pixels
    protect: tuple = ()                # (x, y, r) discs, raw pixels
    paint_mode: str = "fill"           # painted areas: "fill" (inpaint) or "flatten" (to base)
    # Compact islands, markers, texture
    remove_compact: float = None       # mm
    compact_fraction: float = 0.0
    compact_cap: float = 0.15
    markers: tuple = ()                # (x, y) raw pixels
    marker_size_mm: float = 3.0
    marker_height: float = 0.3
    texture_spacing_mm: float = 4.0
    texture_height_mm: float = 0.6
    # Shapes and labels on the image. Each annotation is a dict with raw-pixel
    # coordinates (row 0 at the bottom):
    #   {"type": "label", "x", "y", "text"}          (lines separated by newlines)
    #   {"type": "line", "points": [[x, y], ...], "arrow": bool, "closed": bool}
    #   {"type": "circle", "x", "y", "radius_mm", "filled": bool}
    # Optional per-annotation "width_mm", "height_mm" (shapes) and "mode", "pad"
    # (labels) override the defaults below.
    annotations: tuple = ()
    shape_width_mm: float = 1.2
    shape_height_mm: float = 0.8
    image_label_mode: str = "braille"
    image_label_pad: str = "base"      # "base": flattened to the base; "raised": a flat platform
    image_label_margin_mm: float = 2.0
    # Border and label
    border_width_mm: float = 0.0
    border_height_mm: float = 0.0
    label: tuple = ()
    label_mode: str = "text"
    label_depth_mm: float = None
    label_size_mm: float = 5.0
    label_area_mm: float = 0.0
    label_spacing_mm: float = None

    def effective_clip(self):
        if self.clip is not None:
            return self.clip
        return 0.0 if (self.clean or self.stretch == "tactile") else 1.0

    def to_dict(self):
        return {k: (list(v) if isinstance(v, tuple) else v) for k, v in asdict(self).items()}


@dataclass
class Layers:
    """Hand-painted layers on the mask grid (any may be None) and a version
    number per layer, used to know when cached results are stale."""
    paint: np.ndarray = None           # bool: paint out
    protect: np.ndarray = None         # bool: never mask
    texture: np.ndarray = None         # uint8: 0 none, 1..3 = shaping.TEXTURES
    versions: dict = field(default_factory=dict)


@dataclass
class HeightResult:
    heights: np.ndarray                # final normalised heights, row 0 at the bottom
    factor: int                        # raw pixels per height-map pixel
    pixel_mm: float
    clean: object = None               # CleanResult, if cleaned
    display: np.ndarray = None         # the input at height-map resolution, 0-1, for previews
    image_labels: list = field(default_factory=list)   # (mesh, x, y, level, pad box) per image label
    warnings: list = field(default_factory=list)       # things worth checking before printing
    curve: dict = None                 # brightness rank vs height after the stretch (previews only)
    removed: dict = None               # what source removal took out, for previews (see height_map)


# --- Mask grid --------------------------------------------------------------------

def mask_scale(raw_shape):
    """Raw pixels per mask-grid cell (along each axis)."""
    return max(1, int(np.ceil(max(raw_shape) / MAX_MASK_GRID)))


def mask_grid_shape(raw_shape):
    ms = mask_scale(raw_shape)
    return int(np.ceil(raw_shape[0] / ms)), int(np.ceil(raw_shape[1] / ms))


def regrid(arr, src_scale, dst_scale, out_shape):
    """
    Nearest-neighbour resampling between grids that share an origin: each
    cell of the output (dst_scale raw pixels across) takes the source cell
    (src_scale raw pixels across) under its centre. Scales may be numbers or
    (y, x) pairs.
    """
    sy, sx = (src_scale, src_scale) if np.isscalar(src_scale) else src_scale
    dy, dx = (dst_scale, dst_scale) if np.isscalar(dst_scale) else dst_scale
    iy = np.clip(((np.arange(out_shape[0]) + 0.5) * dy / sy).astype(int), 0, arr.shape[0] - 1)
    ix = np.clip(((np.arange(out_shape[1]) + 0.5) * dx / sx).astype(int), 0, arr.shape[1] - 1)
    return np.asarray(arr)[np.ix_(iy, ix)]


def load_layer(path, raw_shape, integer=False):
    """
    Loads a painted layer saved by the GUI (or any FITS/PNG mask) onto the
    mask grid of an image of raw_shape. Masks of any size are stretched to
    cover the whole image. Non-zero is "on"; with integer=True the values are
    kept (texture maps).
    """
    data = np.nan_to_num(core.load_image(path))
    grid = mask_grid_shape(raw_shape)
    ms = mask_scale(raw_shape)
    src = (raw_shape[0] / data.shape[0], raw_shape[1] / data.shape[1])
    out = regrid(data, src, ms, grid)
    return np.rint(out).astype(np.uint8) if integer else out != 0


def _key(layers, name):
    a = getattr(layers, name) if layers is not None else None
    if a is None or not a.any():
        return None
    return layers.versions.get(name, id(a))


# --- Builder --------------------------------------------------------------------------

class ModelBuilder:
    """
    Builds height maps and meshes for one input image, caching the slow
    stages (star removal and inpainting, noise cleaning) between calls, so an
    interactive preview only redoes what changed. Thread-safe.
    """

    def __init__(self, data, header=None):
        self.raw = np.asarray(data, dtype=np.float32)
        self.header = header
        self.mask_scale = mask_scale(self.raw.shape)
        self.mask_shape = mask_grid_shape(self.raw.shape)
        self._lock = threading.Lock()
        self._cache = {}

    @classmethod
    def from_file(cls, path, hdu=None):
        data, header = core.load_image(path, hdu, return_header=True)
        return cls(data, header)

    def _stage(self, name, key, compute):
        hit = self._cache.get(name)
        if hit is not None and hit[0] == key:
            return hit[1]
        value = compute()
        self._cache[name] = (key, value)
        return value

    def _to_raw(self, layer):
        return regrid(layer, self.mask_scale, 1, self.raw.shape)

    def _sources(self, s, layers):
        paint_key = _key(layers, "paint") if s.paint_mode == "fill" else None
        protect_key = _key(layers, "protect")
        needed = s.mask_stars or s.remove_points or s.masks or paint_key is not None
        if not needed:
            return self.raw, None, None
        key = (s.mask_stars, s.gaia_cache, s.remove_points, s.point_box, s.point_sigma, s.max_point_area,
               tuple(map(tuple, s.masks)), tuple(map(tuple, s.protect)), paint_key, protect_key)

        def compute():
            out, _, parts = sources.remove_sources(
                self.raw, self.header, gaia=s.mask_stars, gaia_cache=s.gaia_cache, points=s.remove_points,
                point_area=s.max_point_area, point_sigma=s.point_sigma, point_box=s.point_box,
                extra_masks=s.masks, protect=s.protect,
                extra_mask=self._to_raw(layers.paint) if paint_key is not None else None,
                protect_mask=self._to_raw(layers.protect) if protect_key is not None else None, parts=True)
            return out, parts
        out, parts = self._stage("sources", key, compute)
        return out, key, parts

    def _cleaned(self, data, src_key, s):
        if not s.clean:
            return data, None
        beam = tuple(s.beam) if s.beam else beam_from_header(self.header)
        key = (src_key, beam, s.detect_sigma, s.peak_sigma, s.min_island,
               tuple(s.smooth_scales) if s.smooth_scales else None, s.noise_box)

        def compute():
            if beam is None:
                logger.warning("No beam given or in the header; assuming 2 px.")
            return clean_image(data, beam_px=beam, smooth_scales_px=s.smooth_scales,
                               detect_sigma=s.detect_sigma, peak_sigma=s.peak_sigma,
                               min_area_px=s.min_island, noise_box_px=s.noise_box)
        result = self._stage("clean", key, compute)
        return result.image, result

    def height_map(self, s, layers=None, with_display=False):
        """The final height map for these settings and painted layers."""
        if s.stretch not in STRETCHES:
            raise ValueError(f"stretch must be one of {', '.join(STRETCHES)}")
        if s.paint_mode not in PAINT_MODES:
            raise ValueError(f"paint_mode must be one of {', '.join(PAINT_MODES)}")
        with self._lock:
            data, src_key, parts = self._sources(s, layers)
            data, clean_result = self._cleaned(data, src_key, s)

        factor = max(int(s.downsample), shaping.auto_downsample_factor(data.shape, s.longest_side_mm, s.min_pixel_mm))
        edge = s.edge if s.edge is not None else (s.detect_sigma if s.clean else None)
        h = core.preprocess_image(
            data, downsample_factor=factor, clip_percentile=s.effective_clip(), nan_value=s.nan_value,
            log_scale=s.stretch == "log", asinh_scale=s.stretch == "asinh", invert=s.invert,
            soften=s.soften, downsample_method=s.downsample_method, ridge=s.ridge,
            tactile=s.stretch == "tactile", tactile_edge=edge, tactile_edge_step=s.edge_step,
            tactile_contrast=s.contrast, tactile_peak_percentile=s.peak_percentile)
        pixel_mm = s.longest_side_mm / max(h.shape) if s.longest_side_mm and s.longest_side_mm > 0 else 1.0
        curve = None
        if with_display:
            src = shaping.block_reduce(data, factor) if factor > 1 else np.asarray(data, np.float32)
            if src.shape == h.shape:
                curve = stretch_curve(src, h)

        if s.remove_compact:
            h, n = shaping.suppress_compact(h, pixel_mm, s.remove_compact,
                                            fraction=s.compact_fraction, cap=s.compact_cap)
            logger.info(f"{'Lowered' if s.compact_fraction > 0 else 'Removed'} {n} compact island(s) "
                        f"under {s.remove_compact:g} mm.")

        def on_grid(name):
            if _key(layers, name) is None:
                return None
            return regrid(getattr(layers, name), self.mask_scale, factor, h.shape)

        if s.paint_mode == "flatten":
            flat = on_grid("paint")
            if flat is not None:
                h = np.where(flat, 0.0, h).astype(np.float32)

        if s.smooth and s.smooth > 0:
            h = ndimage.gaussian_filter(h, s.smooth)

        if s.markers:
            pos = [((x + 0.5) / factor - 0.5, (y + 0.5) / factor - 0.5) for x, y in s.markers]
            h = shaping.add_markers(h, pos, pixel_mm, diameter_mm=s.marker_size_mm, height=s.marker_height)

        tex = on_grid("texture")
        if tex is not None:
            h = shaping.apply_texture(h, tex, pixel_mm, s.texture_spacing_mm,
                                      s.texture_height_mm / s.max_height_mm)

        removed = None
        if with_display and parts is not None:
            # Each removed area on the height-map grid (a model pixel counts if
            # any raw pixel in it was removed), for outlining in previews.
            def to_grid(m):
                if not m.any():
                    return np.zeros(h.shape, bool)
                g = shaping.block_reduce(m.astype(np.float32), factor) > 0 if factor > 1 else m
                return g if g.shape == h.shape else regrid(g, 1, (g.shape[0] / h.shape[0], g.shape[1] / h.shape[1]),
                                                          h.shape)
            removed = {k: to_grid(parts[k]) for k in ("stars", "points", "discs")}
            removed.update(n_stars=parts["n_stars"], n_points=parts["n_points"],
                           fraction=float((parts["stars"] | parts["points"] | parts["discs"]).mean()))

        warnings = []
        h, image_labels = self._annotate(s, h, factor, pixel_mm, warnings)

        display = None
        if with_display:
            d = shaping.block_reduce(self.raw, factor) if factor > 1 else np.array(self.raw, np.float32)
            d = regrid(np.nan_to_num(d, nan=np.nanmin(d)), 1, (d.shape[0] / h.shape[0], d.shape[1] / h.shape[1]),
                       h.shape)
            lo, hi = np.percentile(d, [1, 99.7])
            display = np.arcsinh(np.clip((d - lo) / max(hi - lo, 1e-12), 0, 1) * 10) / np.arcsinh(10)
        return HeightResult(heights=h.astype(np.float32), factor=factor, pixel_mm=pixel_mm,
                            clean=clean_result, display=display, removed=removed, image_labels=image_labels,
                            warnings=warnings, curve=curve)

    PAD_SKIRT_MM = 1.5      # slope from a label pad back up (or down) to the relief around it

    def _annotate(self, s, h, factor, pixel_mm, warnings=None):
        """Raised shapes, then flat pads for the labels on the image (so a leader
        line stops at its label). Returns the heights and the placed labels."""
        if not s.annotations:
            return h, []

        def grid(x, y):
            return (x + 0.5) / factor - 0.5, (y + 0.5) / factor - 0.5

        shapes, labels = [], []
        for a in s.annotations:
            kind = a.get("type")
            common = dict(width_mm=a.get("width_mm", s.shape_width_mm),
                          height=a.get("height_mm", s.shape_height_mm) / s.max_height_mm)
            if kind == "line":
                shapes.append(dict(type="line", points=[grid(*p) for p in a["points"]],
                                   arrow=bool(a.get("arrow")), closed=bool(a.get("closed")), **common))
            elif kind == "circle":
                gx, gy = grid(a["x"], a["y"])
                shapes.append(dict(type="circle", x=gx, y=gy, radius_mm=a["radius_mm"],
                                   filled=bool(a.get("filled")), **common))
            elif kind == "label":
                labels.append(a)
            else:
                raise ValueError(f"Unknown annotation type {kind!r}")
        if shapes:
            narrow = min(sh["width_mm"] for sh in shapes) / pixel_mm
            if narrow < 3 and warnings is not None:
                msg = (f"Raised lines are only {narrow:.1f} model pixels wide ({pixel_mm:.2f} mm pixels), so they "
                       "will be jagged; make them wider, or lower the smallest model pixel / downsampling.")
                warnings.append(msg)
                logger.warning(msg)
            h = shaping.add_shapes(h, shapes, pixel_mm)

        placed = []
        ny, nx = h.shape
        for a in labels:
            mode = a.get("mode", s.image_label_mode)
            pad = a.get("pad", s.image_label_pad)
            if pad not in LABEL_PADS:
                raise ValueError(f"label pad must be one of {', '.join(LABEL_PADS)}")
            lines = [ln for ln in str(a.get("text", "")).split("\n") if ln.strip()]
            if not lines:
                continue
            m, _ = core.generate_label_mesh(lines, mode=mode, depth_mm=s.label_depth_mm if mode == s.label_mode else None,
                                            text_height_mm=s.label_size_mm, line_spacing_mm=s.label_spacing_mm)
            b = m.bounds
            margin = s.image_label_margin_mm
            gx, gy = grid(a["x"], a["y"])
            half_w = (max(abs(b[0]), abs(b[1])) + margin) / pixel_mm
            half_h = (max(abs(b[2]), abs(b[3])) + margin) / pixel_mm
            c0, c1 = int(np.floor(gx - half_w)), int(np.ceil(gx + half_w)) + 1
            r0, r1 = int(np.floor(gy - half_h)), int(np.ceil(gy + half_h)) + 1
            if c0 < 0 or r0 < 0 or c1 > nx or r1 > ny:
                raise ValueError(f"Label {' / '.join(lines)!r} ({(b[1] - b[0]):.0f} x {(b[3] - b[2]):.0f} mm, "
                                 "on its pad) runs off the image; move it inwards"
                                 + (" (braille cannot be scaled down)." if mode == "braille" else "."))
            under = h[r0:r1, c0:c1]
            level = float(under.max()) if pad == "raised" else 0.0
            covered = float(np.mean(np.abs(under - level) > 0.05))
            if covered > 0.05 and warnings is not None:
                msg = (f"Label {' / '.join(lines)!r}: its pad covers relief over {100 * covered:.0f}% of its area; "
                       "move it to clear sky so it doesn't cut into a source.")
                warnings.append(msg)
                logger.warning(msg)
            # Flat pad, with a sloped skirt back to the relief instead of a sheer cut.
            k = max(1, int(np.ceil(self.PAD_SKIRT_MM / pixel_mm)))
            R0, R1, C0, C1 = max(0, r0 - k), min(ny, r1 + k), max(0, c0 - k), min(nx, c1 + k)
            yy, xx = np.mgrid[R0:R1, C0:C1]
            dy = np.maximum(np.maximum(r0 - yy, yy - (r1 - 1)), 0)
            dx = np.maximum(np.maximum(c0 - xx, xx - (c1 - 1)), 0)
            t = np.clip(np.hypot(dx, dy) / (k + 1), 0, 1)
            wgt = t * t * (3 - 2 * t)                     # 0 on the pad, easing to 1 at the skirt edge
            h = h.copy()
            region = h[R0:R1, C0:C1]
            h[R0:R1, C0:C1] = level + (region - level) * wgt
            placed.append((m, gx, gy, level, (c0, r0, c1, r1), " / ".join(lines)))
        return h, placed

    def mesh(self, s, result):
        """The STL mesh (numpy-stl) for a height map from height_map()."""
        return core.generate_mesh(
            result.heights, longest_side_mm=s.longest_side_mm, max_height_mm=s.max_height_mm,
            base_thickness_mm=s.base_thickness_mm, smoothing_sigma=0,
            border_width_mm=s.border_width_mm, border_height_mm=s.border_height_mm,
            label_text=list(s.label), label_mode=s.label_mode, label_depth_mm=s.label_depth_mm,
            label_font_size_mm=s.label_size_mm, label_area_mm=s.label_area_mm,
            label_line_spacing_mm=s.label_spacing_mm,
            image_labels=[(m, x, y, s.base_thickness_mm + level * s.max_height_mm)
                          for m, x, y, level, _, _ in result.image_labels])


def plaque_size(s, result):
    """Approximate overall plaque size (width, depth, height) in mm."""
    ny, nx = result.heights.shape
    p = result.pixel_mm
    bp = int(round(s.border_width_mm / p)) if s.border_width_mm > 0 else 0
    lp = int(round(s.label_area_mm / p)) if s.label_area_mm > 0 else 0
    width = (nx + 2 * bp - 1) * p
    depth = (ny + 2 * bp + lp - 1) * p
    height = s.base_thickness_mm + max(float(result.heights.max()), 0) * s.max_height_mm
    if s.label:
        label_h = s.label_depth_mm if s.label_depth_mm is not None else (
            core.DEFAULT_BRAILLE_DOT_HEIGHT_MM if s.label_mode == "braille" else core.DEFAULT_TEXT_DEPTH_MM)
        height = max(height, s.base_thickness_mm + s.border_height_mm + label_h)
    for m, _, _, level, _, _ in result.image_labels:
        height = max(height, s.base_thickness_mm + level * s.max_height_mm + m.bounds[5])
    return width, depth, height


def fits_bed(size, bed):
    """True if a (width, depth) plaque fits a (width, depth) bed, either way round."""
    w, d = size[:2]
    bw, bd = bed
    return (w <= bw and d <= bd) or (w <= bd and d <= bw)


# --- Templates ---------------------------------------------------------------------------------
#
# Two independent kinds, so any plaque format can be combined with any image recipe:
# a plaque format sets the physical size, border and label strip; an image recipe
# sets how the image is cleaned and stretched. Applying one resets its fields to the
# defaults first, so switching recipes never leaves an option from the previous one on.

PLAQUE_FIELDS = ("longest_side_mm", "max_height_mm", "base_thickness_mm", "border_width_mm",
                 "border_height_mm", "label_area_mm", "label_mode", "label_depth_mm", "label_size_mm",
                 "label_spacing_mm")
RECIPE_FIELDS = ("clip", "smooth", "downsample", "downsample_method", "ridge", "min_pixel_mm", "stretch",
                 "soften", "invert", "edge", "edge_step", "contrast", "peak_percentile", "clean", "beam",
                 "detect_sigma", "peak_sigma", "min_island", "smooth_scales", "noise_box", "mask_stars",
                 "remove_points", "point_box", "point_sigma", "max_point_area", "remove_compact",
                 "compact_fraction", "compact_cap")

PLAQUE_FORMATS = {
    "large": dict(
        short="Large 200 mm",
        title="Large outreach plaque (200 mm)",
        about="For handling at outreach events: 12 mm of relief, a raised border and an 18 mm braille strip "
              "for one line. Fits a 256 x 256 mm bed.",
        settings=dict(longest_side_mm=200, max_height_mm=12, base_thickness_mm=3, border_width_mm=5,
                      border_height_mm=2, label_area_mm=18, label_mode="braille")),
    "standard": dict(
        short="Standard 110 mm",
        title="Standard plaque, two-line label (110 mm)",
        about="The format of the galaxy, Centaurus A and Helix sets: 8 mm of relief and a 24 mm strip for two "
              "braille lines of up to 17 cells.",
        settings=dict(longest_side_mm=110, max_height_mm=8, base_thickness_mm=2, border_width_mm=3,
                      border_height_mm=1.5, label_area_mm=24, label_mode="braille")),
    "small": dict(
        short="Small 90 mm",
        title="Small gift or test print (90 mm)",
        about="Quick to print: 5 mm of relief and a 12 mm strip. Braille can't shrink, so only about 13 cells fit.",
        settings=dict(longest_side_mm=90, max_height_mm=5, base_thickness_mm=2, border_width_mm=3,
                      border_height_mm=1.5, label_area_mm=12, label_mode="braille")),
}

RECIPES = {
    "plain": dict(
        short="Plain",
        title="Plain image (no processing)",
        about="A linear stretch with light smoothing: a starting point for already-clean images.",
        settings=dict()),
    "radio": dict(
        short="Radio map",
        title="Radio map",
        about="Interferometer images (MeerKAT, VLA, ASKAP): noise is masked and empty sky made flat, then an "
              "asinh stretch gives faint lobes real height. The beam is read from the header.",
        settings=dict(clean=True, stretch="asinh")),
    "radio_filaments": dict(
        short="Radio mosaic",
        title="Radio mosaic with thin filaments",
        about="Large mosaics with varying noise and narrow features (e.g. the Galactic Centre): local noise, "
              "a broad detection scale for faint diffuse emission, ridge-preserving downsampling so filaments "
              "stay continuous, and small isolated sources removed.",
        settings=dict(clean=True, stretch="asinh", noise_box=110, smooth_scales=(3.0, 9.0), ridge=0.4,
                      remove_compact=1.3)),
    "galaxy_smooth": dict(
        short="Smooth galaxy",
        title="Optical galaxy, smooth (elliptical, lenticular)",
        about="Foreground stars (Gaia) and compact sources (globular clusters, background galaxies) removed, "
              "so a smooth galaxy feels smooth; then the tactile stretch. Protect the nucleus with the Protect brush.",
        settings=dict(mask_stars=True, remove_points=True, max_point_area=250, stretch="tactile", smooth=1.5)),
    "galaxy_knotty": dict(
        short="Spiral/irregular",
        title="Optical galaxy, spiral or irregular",
        about="Foreground stars removed but compact knots kept, since star-forming knots define these classes; "
              "the tactile stretch with a stronger contrast boost lifts the arms.",
        settings=dict(mask_stars=True, stretch="tactile", contrast=3.0)),
    "nebula": dict(
        short="Nebula",
        title="Nebula with a central star",
        about="Stars removed, the tactile stretch, and small leftovers cleared. Add a Marker at the central "
              "star if it is too faint to feel, and say that it is a marker.",
        settings=dict(mask_stars=True, remove_points=True, stretch="tactile", smooth=2.0, remove_compact=2.0)),
    "press_picture": dict(
        short="Press picture",
        title="Press picture (JWST, Hubble PNG/JPEG)",
        about="Already stretched for display, so heights follow displayed brightness: wide stars removed against "
              "a broad median filter, a strong contrast boost for shells and filaments, and leftover pips cleared.",
        settings=dict(remove_points=True, point_box=21, point_sigma=4.0, max_point_area=800, stretch="tactile",
                      contrast=4.0, remove_compact=3.0)),
}


def apply_template(settings, plaque=None, recipe=None, custom=None):
    """
    A copy of `settings` with a plaque format and/or an image recipe applied
    (names from PLAQUE_FORMATS / RECIPES), and/or a custom template (a dict of
    Settings fields, e.g. loaded from a saved template file).
    """
    d = settings.to_dict()
    base = Settings().to_dict()
    for name, table, fields in ((plaque, PLAQUE_FORMATS, PLAQUE_FIELDS), (recipe, RECIPES, RECIPE_FIELDS)):
        if name is None:
            continue
        if name not in table:
            raise ValueError(f"Unknown template {name!r}; choose from {', '.join(table)}.")
        for f in fields:
            d[f] = base[f]
        d.update(table[name]["settings"])
    if custom:
        unknown = set(custom) - set(base)
        if unknown:
            raise ValueError(f"Unknown settings in template: {', '.join(sorted(unknown))}")
        d.update({k: v for k, v in custom.items() if k not in TEMPLATE_EXCLUDED})
    return Settings(**{k: (tuple(v) if isinstance(v, list) else v) for k, v in d.items()})


# Per-image things a saved template never carries.
TEMPLATE_EXCLUDED = ("annotations", "markers", "masks", "protect", "gaia_cache", "label")


def template_dict(settings):
    """The settings worth saving as a reusable template (no per-image content)."""
    return {k: v for k, v in settings.to_dict().items() if k not in TEMPLATE_EXCLUDED}


# --- Explaining the stretch -----------------------------------------------------------------------

def stretch_curve(stretch_input, heights, n=3000, seed=0):
    """
    A sample of (brightness rank in %, height 0-1) pairs from a preview, for
    charting what the stretch does: how much of the image is flat base, and
    how the brightness range is spread over the relief.
    """
    x = np.asarray(stretch_input, dtype=np.float64).ravel()
    y = np.asarray(heights, dtype=np.float64).ravel()
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    if x.size == 0:
        return None
    rank = np.empty(x.size)
    rank[np.argsort(x, kind="stable")] = np.linspace(0, 100, x.size)
    idx = np.random.default_rng(seed).choice(x.size, size=min(n, x.size), replace=False)
    flat = float(np.mean(y <= 1e-6))
    return dict(rank=rank[idx], height=y[idx], flat_fraction=flat,
                top_fraction=float(np.mean(y >= 0.98)))
