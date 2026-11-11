"""
Shaping the height map for touch: resampling, the tactile stretch, compact
islands and marker domes. Ported from the bvi-models pipelines (MeerKAT
Galactic Centre, galaxy morphology, Helix Nebula).
"""
import logging

import numpy as np
from scipy import ndimage

logger = logging.getLogger(__name__)


# --- Resampling ---------------------------------------------------------------------

def block_reduce(data, factor, ridge=0.0):
    """
    Downsamples by an integer factor with a block mean: a proper low-pass
    filter, so noise averages away instead of aliasing into grit (as linear
    interpolation does).

    ridge mixes in that fraction of the block maximum: (1 - ridge) * mean +
    ridge * max. This keeps thin features (radio filaments, narrow arms) as
    continuous ridges a finger can follow instead of diluting them away. Use
    it only on maps whose noise is already masked, or the maximum adds noise.

    NaN-aware: a block less than half valid becomes NaN. Rows and columns that
    do not fill a whole block are cropped from the top and right edges.
    """
    k = int(factor)
    if k <= 1:
        return np.array(data, dtype=np.float32, copy=True)
    a = np.asarray(data, dtype=np.float64)
    ny, nx = (a.shape[0] // k) * k, (a.shape[1] // k) * k
    if ny == 0 or nx == 0:
        raise ValueError(f"Downsampling by {k} leaves an empty image ({a.shape[1]}x{a.shape[0]} px).")
    if (ny, nx) != a.shape:
        logger.info(f"Block downsampling crops {a.shape[0] - ny} row(s) and {a.shape[1] - nx} column(s).")
    b = a[:ny, :nx].reshape(ny // k, k, nx // k, k)
    valid = np.isfinite(b)
    count = valid.sum(axis=(1, 3))
    with np.errstate(invalid="ignore", divide="ignore"):
        out = np.where(valid, b, 0.0).sum(axis=(1, 3)) / count
    if ridge > 0:
        mx = np.where(valid, b, -np.inf).max(axis=(1, 3))
        out = (1 - ridge) * out + ridge * mx
    out[count < 0.5 * k * k] = np.nan
    return out.astype(np.float32)


def auto_downsample_factor(shape, longest_side_mm, min_pixel_mm):
    """
    The integer factor that makes model pixels at least min_pixel_mm across.
    Pixels much finer than a 0.4 mm nozzle can print only bloat the STL (and
    slow the slicer) without adding anything a finger can feel.
    """
    if not longest_side_mm or longest_side_mm <= 0 or not min_pixel_mm or min_pixel_mm <= 0:
        return 1
    return max(1, int(np.floor(min_pixel_mm * max(shape) / longest_side_mm)))


# --- Tactile stretch ------------------------------------------------------------------

def tactile_stretch(img, edge, peak_percentile=99.9, edge_step=0.2, contrast=1.5, contrast_scale=0.025):
    """
    Log brightness stretch designed for touch (morph.to_height in bvi-models).

    - Everything fainter than `edge` (the outer isophote, in data units) is
      exactly 0, so it prints as flat base.
    - The edge itself maps to `edge_step` of full height, so the outline of
      the object is a definite step a finger finds at once.
    - The peak is capped at the `peak_percentile` brightness, so a bright core
      becomes a small rounded summit instead of taking most of the height
      range and leaving discs and arms flat.
    - A local-contrast boost in log space (strength `contrast`, on structures
      up to `contrast_scale` of the map size) lifts arms, bars, knots and
      shells above their surroundings; smooth profiles barely change.

    Returns heights in [0, 1].
    """
    img = np.asarray(img, dtype=np.float64)
    inside = img > edge
    if not inside.any():
        raise ValueError(f"Nothing in the image is brighter than the tactile edge ({edge:.4g}). "
                         "Lower --edge.")
    peak = np.percentile(img[inside], peak_percentile)
    if not peak > edge:
        raise ValueError(f"The image peak ({peak:.4g}) is not above the tactile edge ({edge:.4g}).")
    img = np.minimum(img, peak)
    lg = np.log10(np.clip(img, edge, None) / edge) / np.log10(peak / edge)
    if contrast:
        sigma = contrast_scale * max(img.shape)
        # Normalised blur, so the object's edge does not pull the boost down.
        w = ndimage.gaussian_filter(inside.astype(float), sigma)
        local = ndimage.gaussian_filter(np.where(inside, lg, 0.0), sigma) / np.maximum(w, 1e-6)
        lg = np.clip(lg + contrast * (lg - local), 0, None)
    h = np.where(inside, edge_step + (1 - edge_step) * lg, 0.0)
    return np.clip(h / h.max(), 0, 1).astype(np.float32)


# --- Compact islands --------------------------------------------------------------------

def suppress_compact(heights, pixel_mm, min_size_mm, fraction=0.0, cap=0.0):
    """
    Removes (or lowers) isolated islands smaller than min_size_mm across on the
    model: background radio galaxies, compact HII regions and the like, which
    print as 1-2 mm pips that compete with the story and can be mistaken for
    braille dots.

    Islands are connected areas above zero, so this needs a map with a flat
    zero background (after --clean or --tactile). With fraction > 0 they are
    kept as low pips instead: lowered to `fraction` of their height and
    capped at `cap` of full height.

    Returns (heights, number of islands affected).
    """
    h = np.array(heights, dtype=np.float32, copy=True)
    on = h > 0
    if on.mean() > 0.95:
        logger.warning("Compact-island removal needs a flat, zero background (use --clean or "
                       "--tactile); almost the whole map is above zero, so nothing was removed.")
        return h, 0
    labels, n = ndimage.label(on)
    if n == 0:
        return h, 0
    area = ndimage.sum(on, labels, np.arange(1, n + 1))
    min_area = (min_size_mm / pixel_mm) ** 2
    small = np.nonzero(area < min_area)[0] + 1
    compact = np.isin(labels, small)
    h[compact] = np.minimum(h[compact] * fraction, cap)
    return h, int(small.size)


# --- Markers -------------------------------------------------------------------------------

def add_markers(heights, positions_px, pixel_mm, diameter_mm=3.0, height=0.3):
    """
    Adds a smooth round dome at each (x, y) pixel position, standing `height`
    (a fraction of full relief) above the local relief around it. For things
    too faint to feel but important to the story, such as a nebula's central
    star. A marker is not data: say so when presenting the model.
    """
    h = np.array(heights, dtype=np.float32, copy=True)
    r_px = diameter_mm / 2 / pixel_mm
    if r_px < 1:
        logger.warning(f"Marker ({diameter_mm:g} mm) is smaller than one model pixel ({pixel_mm:.2f} mm).")
    yy, xx = np.indices(h.shape)
    for x, y in positions_px:
        if not (0 <= x < h.shape[1] and 0 <= y < h.shape[0]):
            raise ValueError(f"Marker at pixel ({x:.0f}, {y:.0f}) is outside the image.")
        rr = np.hypot(yy - y, xx - x) / max(r_px, 1e-6)
        dome = np.sqrt(np.clip(1 - rr ** 2, 0, 1))
        ring = (rr > 1.5) & (rr < 3)
        base = float(np.median(h[ring])) if ring.any() else float(h[int(round(y)), int(round(x))])
        h = np.where(rr < 1, np.maximum(h, base + height * dome), h)
    return np.clip(h, 0, 1)


# --- Tactile textures ---------------------------------------------------------------------

# Values in a texture map: 0 = none, then these patterns in order.
TEXTURES = ("dots", "lines", "grid")


def _ridges(x_mm, y_mm, spacing_mm, angle_deg):
    """Parallel rounded ridges, 0.4 x spacing wide, across the given angle."""
    a = np.radians(angle_deg)
    u = -x_mm * np.sin(a) + y_mm * np.cos(a)
    d = np.abs(u - np.round(u / spacing_mm) * spacing_mm)
    half = 0.2 * spacing_mm
    return np.sqrt(np.clip(1 - (d / half) ** 2, 0, 1))


def texture_pattern(shape, pixel_mm, kind, spacing_mm):
    """
    A tactile fill pattern in [0, 1] over the whole grid, anchored to the
    model (so neighbouring regions with the same pattern line up):

    - "dots": round bumps 0.5 x spacing across on a hexagonal lattice. Keep
      the spacing well away from braille (2.5 mm dots, 1.5 mm across) so the
      texture is never mistaken for writing.
    - "lines": parallel rounded ridges at 45 degrees.
    - "grid": ridges in both directions (a crosshatch).
    """
    yy, xx = np.indices(shape, dtype=np.float64) * pixel_mm
    s = float(spacing_mm)
    if kind == "dots":
        row_h = s * np.sqrt(3) / 2
        r_idx = np.round(yy / row_h)
        cy = r_idx * row_h
        off = (r_idx % 2) * s / 2
        cx = np.round((xx - off) / s) * s + off
        R = 0.25 * s
        return np.sqrt(np.clip(1 - (np.hypot(xx - cx, yy - cy) / R) ** 2, 0, 1))
    if kind == "lines":
        return _ridges(xx, yy, s, 45.0)
    if kind == "grid":
        return np.maximum(_ridges(xx, yy, s, 0.0), _ridges(xx, yy, s, 90.0))
    raise ValueError(f"Unknown texture {kind!r}; choose from {', '.join(TEXTURES)}.")


def apply_texture(heights, texture_map, pixel_mm, spacing_mm, height):
    """
    Raises a tactile pattern on top of the relief wherever texture_map is
    non-zero (1 = dots, 2 = lines, 3 = grid; see TEXTURES). `height` is the
    pattern's height as a fraction of full relief. Textured summits may go
    slightly above 1, so the pattern is never flattened off the top.

    Tactile-graphics guidance: textures need elements at least ~0.5 mm high
    and 2-5 mm apart, and should differ clearly from each other and from
    braille.
    """
    h = np.array(heights, dtype=np.float32, copy=True)
    tm = np.asarray(texture_map)
    if spacing_mm < 6 * pixel_mm:
        logger.warning(f"Texture spacing {spacing_mm:g} mm is only {spacing_mm / pixel_mm:.1f} model pixels; "
                       "the pattern will be coarse. Use a larger spacing or smaller --min-pixel.")
    for value, kind in enumerate(TEXTURES, start=1):
        m = tm == value
        if m.any():
            h[m] += height * texture_pattern(h.shape, pixel_mm, kind, spacing_mm)[m].astype(np.float32)
    return h


# --- Display ---------------------------------------------------------------------------------

def hillshade(heights, pixel_mm=1.0, max_height_mm=10.0, azimuth_deg=315.0, altitude_deg=45.0):
    """
    Shaded relief of a height map (rows top to bottom as displayed), lit from
    the north-west: what the model will look like, for previews.
    """
    z = np.asarray(heights, dtype=np.float64) * max_height_mm
    gy, gx = np.gradient(z, pixel_mm)
    slope = np.pi / 2 - np.arctan(np.hypot(gx, gy))
    aspect = np.arctan2(-gx, gy)
    az = np.radians(360.0 - azimuth_deg + 90.0)
    alt = np.radians(altitude_deg)
    shaded = np.sin(alt) * np.sin(slope) + np.cos(alt) * np.cos(slope) * np.cos(az - aspect)
    return np.clip(shaded, 0, 1)


# --- Raised shapes (leader lines, arrows, rings) ---------------------------------------------

def _segment_distance(xx, yy, p0, p1):
    """Distance from each (xx, yy) to the segment p0-p1."""
    (x0, y0), (x1, y1) = p0, p1
    dx, dy = x1 - x0, y1 - y0
    L2 = dx * dx + dy * dy
    t = np.zeros_like(xx) if L2 == 0 else np.clip(((xx - x0) * dx + (yy - y0) * dy) / L2, 0, 1)
    return np.hypot(xx - (x0 + t * dx), yy - (y0 + t * dy))


def _in_triangle(xx, yy, a, b, c):
    def side(p, q):
        return (xx - q[0]) * (p[1] - q[1]) - (p[0] - q[0]) * (yy - q[1])
    d1, d2, d3 = side(a, b), side(b, c), side(c, a)
    neg = (d1 < 0) | (d2 < 0) | (d3 < 0)
    pos = (d1 > 0) | (d2 > 0) | (d3 > 0)
    return ~(neg & pos)


def add_shapes(heights, shapes, pixel_mm):
    """
    Raises tactile shapes on top of the relief, following its surface (like
    textures). Each shape is a dict in height-map pixel coordinates (x =
    column, y = row, row 0 at the bottom), with "width_mm" and "height"
    (a fraction of full relief):

    - {"type": "line", "points": [(x, y), ...], "closed": bool, "arrow": bool}
      A rounded ridge along the polyline (closed: back to the start), with
      a filled arrowhead at the last point if "arrow".
    - {"type": "circle", "x", "y", "radius_mm", "filled": bool}
      A ring, or a flat-topped disc.

    Where shapes overlap, the higher one wins (heights are not added up).
    Tactile-graphics guidance: raised lines about 1 mm or more wide and
    0.5-1 mm high; arrowheads clearly larger than the line.
    """
    h = np.array(heights, dtype=np.float32, copy=True)
    add = np.zeros_like(h)
    ny, nx = h.shape

    def stamp(bbox, profile_fn, height):
        x0, y0, x1, y1 = bbox
        c0, c1 = max(0, int(np.floor(x0))), min(nx, int(np.ceil(x1)) + 1)
        r0, r1 = max(0, int(np.floor(y0))), min(ny, int(np.ceil(y1)) + 1)
        if c0 >= c1 or r0 >= r1:
            return
        yy, xx = np.mgrid[r0:r1, c0:c1].astype(np.float64)
        prof = profile_fn(xx, yy)
        np.maximum(add[r0:r1, c0:c1], height * prof, out=add[r0:r1, c0:c1])

    for sh in shapes:
        w = sh.get("width_mm", 1.2) / pixel_mm
        half = w / 2
        height = sh.get("height", 0.1)

        def ridge(d, half=half):
            return np.sqrt(np.clip(1 - (d / max(half, 0.5)) ** 2, 0, 1))

        if sh["type"] == "line":
            pts = [tuple(map(float, p)) for p in sh["points"]]
            if len(pts) < 2:
                continue
            if sh.get("closed") and len(pts) > 2:
                pts = pts + [pts[0]]
            head = None
            if sh.get("arrow"):
                (xa, ya), (xb, yb) = pts[-2], pts[-1]
                L = np.hypot(xb - xa, yb - ya)
                if L > 0:
                    ux, uy = (xb - xa) / L, (yb - ya) / L
                    hl = min(3.5 * w, 0.8 * L)           # arrowhead length
                    hw = 1.6 * w                         # arrowhead half-width
                    base = (xb - ux * hl, yb - uy * hl)
                    head = ((xb, yb), (base[0] - uy * hw, base[1] + ux * hw), (base[0] + uy * hw, base[1] - ux * hw))
                    pts[-1] = (base[0] + ux * 0.3 * hl, base[1] + uy * 0.3 * hl)
            for p0, p1 in zip(pts[:-1], pts[1:]):
                bbox = (min(p0[0], p1[0]) - half - 1, min(p0[1], p1[1]) - half - 1,
                        max(p0[0], p1[0]) + half + 1, max(p0[1], p1[1]) + half + 1)
                stamp(bbox, lambda xx, yy, p0=p0, p1=p1: ridge(_segment_distance(xx, yy, p0, p1)), height)
            if head is not None:
                xs, ys = [p[0] for p in head], [p[1] for p in head]
                bbox = (min(xs) - 1, min(ys) - 1, max(xs) + 1, max(ys) + 1)
                stamp(bbox, lambda xx, yy, t=head: _in_triangle(xx + 0.0, yy + 0.0, *t).astype(float), height)
        elif sh["type"] == "circle":
            cx, cy = float(sh["x"]), float(sh["y"])
            R = sh["radius_mm"] / pixel_mm
            bbox = (cx - R - half - 1, cy - R - half - 1, cx + R + half + 1, cy + R + half + 1)
            if sh.get("filled"):
                stamp(bbox, lambda xx, yy: ridge(np.clip(np.hypot(xx - cx, yy - cy) - R, 0, None)), height)
            else:
                stamp(bbox, lambda xx, yy: ridge(np.abs(np.hypot(xx - cx, yy - cy) - R)), height)
        else:
            raise ValueError(f"Unknown shape type {sh['type']!r}")
    return h + add
