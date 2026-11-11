#!/usr/bin/env python3
# AstroTouch: Convert astronomical FITS images to 3D printable STL files.
# Designed primarily for astronomy outreach for the blind and visually impaired.
#
# Original by Enzo Peres Afonso 2025
# Refactored and Enhanced

import argparse
import logging
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits
from stl import mesh
import pyvista as pv
from scipy.ndimage import gaussian_filter, zoom

try:
    from . import shaping
except ImportError:  # imported as a top-level module (src/astrotouch on sys.path)
    import shaping

# Try to import tqdm for progress bars, but don't fail if not present
try:
    from tqdm import tqdm
except ImportError:
    tqdm = None

# Logging is configured in main(), not on import, so that importing this module
# from another script does not reconfigure that script's logging.
logger = logging.getLogger(__name__)

# Standard Grade 1 Braille mapping (a-z, 0-9, space, UEB punctuation)
BRAILLE_PATTERNS = {
    'a': [1,0,0,0,0,0], 'b': [1,1,0,0,0,0], 'c': [1,0,0,1,0,0], 'd': [1,0,0,1,1,0],
    'e': [1,0,0,0,1,0], 'f': [1,1,0,1,0,0], 'g': [1,1,0,1,1,0], 'h': [1,1,0,0,1,0],
    'i': [0,1,0,1,0,0], 'j': [0,1,0,1,1,0], 'k': [1,0,1,0,0,0], 'l': [1,1,1,0,0,0],
    'm': [1,0,1,1,0,0], 'n': [1,0,1,1,1,0], 'o': [1,0,1,0,1,0], 'p': [1,1,1,1,0,0],
    'q': [1,1,1,1,1,0], 'r': [1,1,1,0,1,0], 's': [0,1,1,1,0,0], 't': [0,1,1,1,1,0],
    'u': [1,0,1,0,0,1], 'v': [1,1,1,0,0,1], 'w': [0,1,0,1,1,1], 'x': [1,0,1,1,0,1],
    'y': [1,0,1,1,1,1], 'z': [1,0,1,0,1,1], ' ': [0,0,0,0,0,0],
    '1': [1,0,0,0,0,0], '2': [1,1,0,0,0,0], '3': [1,0,0,1,0,0], '4': [1,0,0,1,1,0],
    '5': [1,0,0,0,1,0], '6': [1,1,0,1,0,0], '7': [1,1,0,1,1,0], '8': [1,1,0,0,1,0],
    '9': [0,1,0,1,0,0], '0': [0,1,0,1,1,0], 
    '#': [0,0,1,1,1,1], # Number sign
    '^': [0,0,0,0,0,1], # Capital sign (Dot 6)
    '.': [0,1,0,0,1,1], ',': [0,1,0,0,0,0], '-': [0,0,1,0,0,1],
    '!': [0,1,1,0,1,0], '?': [0,1,1,0,0,1], ':': [0,1,0,0,1,0], ';': [0,1,1,0,0,0]
}

# UEB grade 1 encoding. Cells are 6-element dot lists [d1, d2, d3, d4, d5, d6].
_LETTERS = {c: BRAILLE_PATTERNS[c] for c in "abcdefghijklmnopqrstuvwxyz"}
_SPACE = [0, 0, 0, 0, 0, 0]
_CAPITAL = [0, 0, 0, 0, 0, 1]        # dot 6
_NUMBER = [0, 0, 1, 1, 1, 1]         # dots 3456
_GRADE1 = [0, 0, 0, 0, 1, 1]         # dots 56: needed before a-j right after a number
_DIGITS = "jabcdefghi"               # UEB digits 0-9 are the letters j, a-i
_PUNCTUATION = {
    '.': [[0, 1, 0, 0, 1, 1]], ',': [[0, 1, 0, 0, 0, 0]], '-': [[0, 0, 1, 0, 0, 1]],
    '!': [[0, 1, 1, 0, 1, 0]], '?': [[0, 1, 1, 0, 0, 1]], ':': [[0, 1, 0, 0, 1, 0]],
    ';': [[0, 1, 1, 0, 0, 0]], "'": [[0, 0, 1, 0, 0, 0]], '\u2019': [[0, 0, 1, 0, 0, 0]],
    '/': [[0, 0, 0, 1, 1, 1], [0, 0, 1, 1, 0, 0]],
    '(': [[0, 0, 0, 0, 1, 0], [1, 1, 0, 0, 0, 1]],
    ')': [[0, 0, 0, 0, 1, 0], [0, 0, 1, 1, 1, 0]],
}

# Default physical label heights (mm). Standard braille dots are ~0.5-0.9 mm high.
DEFAULT_TEXT_DEPTH_MM = 1.5
DEFAULT_BRAILLE_DOT_HEIGHT_MM = 0.8


def _is_letter(c):
    return c.isascii() and c.isalpha()


def ueb_cells(text):
    """
    Encodes text as Unified English Braille grade 1 (uncontracted) cells.

    - A run of two or more letters that are all capitals ("NGC") gets the
      capitals-word indicator; any other capital gets its own capital sign, so
      "SBb" is cap-s cap-b b (the capitals-word indicator would make it "SBB").
    - Digits get the number sign and use the letters a-j; a period or comma
      between digits stays in numeric mode; a letter a-j straight after a
      number gets the grade 1 indicator so it is not read as a digit.

    Raises ValueError for characters that cannot be encoded, rather than
    silently dropping them.
    """
    unsupported = sorted({c for c in text
                          if not (c == ' ' or c in "0123456789" or _is_letter(c) or c in _PUNCTUATION)})
    if unsupported:
        raise ValueError(f"Cannot write {', '.join(repr(c) for c in unsupported)} in braille "
                         f"(label: {text!r}). Remove or replace these characters.")

    cells = []
    in_number = False
    i, n = 0, len(text)
    while i < n:
        c = text[i]
        if c == ' ':
            cells.append(_SPACE)
            in_number = False
            i += 1
        elif c in "0123456789":
            if not in_number:
                cells.append(_NUMBER)
                in_number = True
            cells.append(_LETTERS[_DIGITS[int(c)]])
            i += 1
        elif c in ".," and in_number and i + 1 < n and text[i + 1] in "0123456789":
            cells += _PUNCTUATION[c]      # decimal point / digit-group comma
            i += 1
        elif _is_letter(c):
            j = i
            while j < n and _is_letter(text[j]):
                j += 1
            run = text[i:j]
            if len(run) > 1 and run.isupper():
                cells += [_CAPITAL, _CAPITAL] + [_LETTERS[ch.lower()] for ch in run]
            else:
                for k, ch in enumerate(run):
                    if ch.isupper():
                        cells.append(_CAPITAL)
                    elif k == 0 and in_number and ch in "abcdefghij":
                        cells.append(_GRADE1)
                    cells.append(_LETTERS[ch.lower()])
            in_number = False
            i = j
        else:
            cells += _PUNCTUATION[c]
            in_number = False
            i += 1
    return cells


BRAILLE_DOT_RADIUS_MM = 0.75
BRAILLE_DOT_SPACING_MM = 2.5
BRAILLE_CELL_SPACING_MM = 6.0
BRAILLE_LINE_SPACING_MM = 10.0   # standard braille line pitch
TEXT_LINE_SPACING = 1.7          # text line pitch, as a multiple of the cap height


def generate_braille_mesh(text, dot_radius=BRAILLE_DOT_RADIUS_MM, dot_height=DEFAULT_BRAILLE_DOT_HEIGHT_MM,
                          dot_spacing=BRAILLE_DOT_SPACING_MM, cell_spacing=BRAILLE_CELL_SPACING_MM):
    """Generates a 3D mesh of Braille dots for the given text (UEB grade 1).

    Returns None if the text produces no dots (e.g. only spaces).
    """
    meshes = []

    # 6-dot grid offsets (x, y)
    dot_offsets = [
        (0, 2*dot_spacing), (0, dot_spacing), (0, 0),       # Dots 1, 2, 3
        (dot_spacing, 2*dot_spacing), (dot_spacing, dot_spacing), (dot_spacing, 0) # Dots 4, 5, 6
    ]

    for cell_index, pattern in enumerate(ueb_cells(text)):
        x_pos = cell_index * cell_spacing
        for i, bit in enumerate(pattern):
            if bit:
                dx, dy = dot_offsets[i]
                dot = pv.Sphere(radius=dot_radius, center=(x_pos + dx, dy, 0))
                dot.scale([1, 1, dot_height/dot_radius], inplace=True)
                dot = dot.clip(normal='-z', origin=(0, 0, 0)).fill_holes(1000)
                meshes.append(dot)

    if not meshes:
        return None

    return meshes[0].merge(meshes[1:]) if len(meshes) > 1 else meshes[0]


_CAP_METRICS = None


def _vector_text(text):
    """Raw vtkVectorText outline: baseline near y=0, not centred, flat at z=0."""
    from vtkmodules.vtkRenderingFreeType import vtkVectorText
    src = vtkVectorText()
    src.SetText(text)
    return src


def _cap_metrics():
    """(bottom, height) of a capital H in vtkVectorText units, used as the cap band."""
    global _CAP_METRICS
    if _CAP_METRICS is None:
        src = _vector_text("H")
        src.Update()
        b = src.GetOutput().GetBounds()
        _CAP_METRICS = (b[2], b[3] - b[2])
    return _CAP_METRICS


def generate_text_mesh(text, cap_height_mm=5.0, depth_mm=DEFAULT_TEXT_DEPTH_MM):
    """
    Raised print text: capitals are cap_height_mm tall, the cap band is centred
    on y=0 (descenders hang below it, as in print) and the letters stand from
    z=0 to z=depth_mm. Returns None for blank text.

    (pyvista's Text3D centres the bounding box in x, y and z, which sinks half
    the depth into the surface and moves lines with descenders.)
    """
    if not text.strip():
        return None
    from vtkmodules.vtkFiltersModeling import vtkLinearExtrusionFilter
    from vtkmodules.vtkFiltersCore import vtkTriangleFilter
    src = _vector_text(text)
    extrude = vtkLinearExtrusionFilter()
    extrude.SetInputConnection(src.GetOutputPort())
    extrude.SetExtrusionTypeToVectorExtrusion()
    extrude.SetVector(0, 0, 1)
    extrude.SetScaleFactor(depth_mm)
    extrude.CappingOn()
    tri = vtkTriangleFilter()
    tri.SetInputConnection(extrude.GetOutputPort())
    tri.Update()
    m = pv.wrap(tri.GetOutput()).copy()
    cap_bottom, cap_h = _cap_metrics()
    s = cap_height_mm / cap_h
    pts = np.asarray(m.points, dtype=float)
    pts[:, 0] *= s
    pts[:, 1] = (pts[:, 1] - cap_bottom - cap_h / 2) * s
    m.points = pts
    return m


def _label_lines(label_text):
    """Accept a string (lines split on newlines) or a list of strings."""
    if not label_text:
        return []
    if isinstance(label_text, str):
        label_text = label_text.split("\n")
    return [line.strip() for line in label_text if line and line.strip()]


def generate_label_mesh(lines, mode="text", depth_mm=None, text_height_mm=5.0,
                        line_spacing_mm=None, max_width_mm=None):
    """
    A label of one or more lines as a single mesh, centred on x=0 and y=0 and
    standing on z=0.

    - Braille (UEB grade 1) keeps its standard physical size; a line wider than
      max_width_mm raises ValueError, because shrunken braille is unreadable.
    - Text lines share one cap height (text_height_mm); if any line is wider
      than max_width_mm, all lines shrink by the same factor.
    - line_spacing_mm is the distance between line centres (default 10 mm for
      braille, 1.7 x the cap height for text).

    Returns (mesh, block_height_mm), or (None, 0) if there is nothing to write.
    """
    lines = _label_lines(lines)
    if not lines:
        return None, 0.0
    if mode not in ("text", "braille"):
        raise ValueError(f"label mode must be 'text' or 'braille', not {mode!r}")

    if mode == "braille":
        depth = DEFAULT_BRAILLE_DOT_HEIGHT_MM if depth_mm is None else depth_mm
        pitch = BRAILLE_LINE_SPACING_MM if line_spacing_mm is None else line_spacing_mm
        line_h = 2 * BRAILLE_DOT_SPACING_MM + 2 * BRAILLE_DOT_RADIUS_MM
        parts = []
        for line in lines:
            n_cells = len(ueb_cells(line))
            width = (n_cells - 1) * BRAILLE_CELL_SPACING_MM + BRAILLE_DOT_SPACING_MM + 2 * BRAILLE_DOT_RADIUS_MM
            if max_width_mm is not None and width > max_width_mm:
                raise ValueError(
                    f"Braille label line {line!r} is {width:.0f} mm wide ({n_cells} cells) but only "
                    f"{max_width_mm:.0f} mm fits. Shorten it, split it over more lines, or make the "
                    "model larger: braille cannot be scaled down.")
            m = generate_braille_mesh(line, dot_height=depth)
            if m is None:
                continue
            # Centre on the cell grid (not the dots, which vary by letter):
            # dot columns span x = 0..dot_spacing per cell, rows y = 0..2*dot_spacing.
            m.translate([-((n_cells - 1) * BRAILLE_CELL_SPACING_MM + BRAILLE_DOT_SPACING_MM) / 2,
                         -BRAILLE_DOT_SPACING_MM, 0], inplace=True)
            parts.append(m)
    else:
        depth = DEFAULT_TEXT_DEPTH_MM if depth_mm is None else depth_mm
        cap = text_height_mm
        parts = [generate_text_mesh(line, cap_height_mm=1.0, depth_mm=depth) for line in lines]
        widest = max(p.bounds[1] - p.bounds[0] for p in parts)
        if max_width_mm is not None and widest * cap > max_width_mm:
            new_cap = max_width_mm / widest
            logger.warning(f"Text label too wide at {cap:g} mm cap height; "
                           f"shrinking all lines to {new_cap:.2f} mm to fit.")
            cap = new_cap
        pitch = TEXT_LINE_SPACING * text_height_mm if line_spacing_mm is None else line_spacing_mm
        line_h = 1.3 * cap  # cap band plus descenders
        for p in parts:
            pts = np.asarray(p.points, dtype=float)
            pts[:, :2] *= cap
            pts[:, 0] -= (pts[:, 0].min() + pts[:, 0].max()) / 2
            p.points = pts

    # Stack lines top to bottom, the block centred on y=0.
    n = len(parts)
    for k, p in enumerate(parts):
        p.translate([0, ((n - 1) / 2 - k) * pitch, 0], inplace=True)
    block_h = (n - 1) * pitch + line_h
    merged = parts[0].merge(parts[1:]) if n > 1 else parts[0]
    return merged, block_h


def load_fits_data(filepath, hdu_index=None, return_header=False):
    """
    Loads image data from a FITS file.
    If hdu_index is None, it attempts to find the first 2D image HDU.
    With return_header=True, returns (data, header) of that HDU.
    """
    logger.info(f"Loading FITS file: {filepath}")
    try:
        with fits.open(filepath) as hdul:
            if hdu_index is not None:
                if hdu_index >= len(hdul):
                    raise ValueError(f"HDU index {hdu_index} out of range. Max index: {len(hdul)-1}")
                hdu = hdul[hdu_index]
                if hdu.data is None or hdu.data.ndim != 2:
                    raise ValueError(f"HDU {hdu_index} does not contain 2D image data.")
            else:
                # Auto-detect HDU
                hdu = None
                for i, h in enumerate(hdul):
                    if h.data is not None and h.data.ndim == 2:
                        logger.info(f"Auto-detected 2D image data in HDU {i}")
                        hdu = h
                        break
                if hdu is None:
                    raise ValueError("Could not find any 2D image data in the FITS file.")
            
            # Return copy to avoid issues with closed memmap
            data = hdu.data.astype(np.float32).copy()
            return (data, hdu.header.copy()) if return_header else data
    except Exception as e:
        logger.error(f"Could not read FITS file: {e}")
        raise

IMAGE_SUFFIXES = {".png", ".jpg", ".jpeg", ".tif", ".tiff", ".bmp"}


def load_image(filepath, hdu_index=None, return_header=False):
    """
    Loads a FITS file (see load_fits_data) or an ordinary image (PNG, JPEG,
    TIFF), such as a press-release picture. Colour images are reduced to the
    mean of their channels. Images are flipped so that, as for FITS, row 0 is
    the bottom of the picture; the model then matches the picture. Ordinary
    images have no header (None).
    """
    if Path(filepath).suffix.lower() not in IMAGE_SUFFIXES:
        return load_fits_data(filepath, hdu_index, return_header=return_header)
    from PIL import Image
    logger.info(f"Loading image file: {filepath}")
    with Image.open(filepath) as im:
        if im.mode in ("I;16", "I;16B", "I;16L", "I", "F", "L"):
            data = np.asarray(im, dtype=np.float32)
        else:
            data = np.asarray(im.convert("RGB"), dtype=np.float32).mean(axis=2)
    data = np.ascontiguousarray(np.flipud(data))
    logger.info("Display images were already stretched for the picture, so heights follow the "
                "displayed brightness; a linear stretch or --tactile usually suits them best.")
    return (data, None) if return_header else data


def robust_rms(image_data):
    """
    Median-absolute-deviation noise estimate, ignoring the brightest pixels so a
    strong source does not inflate it.
    """
    finite = image_data[np.isfinite(image_data)]
    if finite.size == 0:
        raise ValueError("Image contains no finite values.")
    quiet = finite[np.abs(finite) < np.percentile(np.abs(finite), 90)]
    if quiet.size == 0:
        quiet = finite
    rms = 1.4826 * np.median(np.abs(quiet - np.median(quiet)))
    return float(rms) if rms > 0 else float(np.std(finite))


def preprocess_image(
    image_data,
    downsample_factor=1,
    clip_percentile=1.0,
    nan_value=None,
    log_scale=False,
    asinh_scale=False,
    invert=False,
    soften=0.0,
    downsample_method="mean",
    ridge=0.0,
    tactile=False,
    tactile_edge=None,
    tactile_edge_step=0.2,
    tactile_contrast=1.5,
    tactile_peak_percentile=99.9
):
    """
    Applies various processing steps to the image data.
    Returns normalized data in range [0, 1]. The input array is not modified.

    downsample_method: "mean" (block mean, a proper low-pass filter; with
        `ridge`, mixes in that fraction of the block maximum to keep thin
        ridges) or "linear" (the old linear interpolation).
    soften: divide by this multiple of the noise rms before --log/--asinh, so
        the stretch compresses whatever the data's units (see below).
    tactile: use shaping.tactile_stretch instead of log/asinh. tactile_edge is
        the outer isophote in data units (default: sky level + 3 x noise).
    """
    if downsample_method not in ("mean", "linear"):
        raise ValueError(f"downsample_method must be 'mean' or 'linear', not {downsample_method!r}")
    if sum(map(bool, (log_scale, asinh_scale, tactile))) > 1:
        raise ValueError("Choose only one of log, asinh and tactile scaling.")
    # Work on a copy: callers (e.g. the GUI) cache the raw data and reuse it.
    image_data = np.array(image_data, dtype=np.float32, copy=True)
    finite_mask = np.isfinite(image_data)
    if not finite_mask.any():
        raise ValueError("Image contains no finite pixels (all NaN/inf).")

    noise = robust_rms(image_data) if soften > 0 else None

    if downsample_factor > 1 and downsample_method == "mean":
        logger.info(f"Downsampling by factor {downsample_factor} (block mean"
                    + (f" + {ridge:g} block max)..." if ridge > 0 else ")..."))
        image_data = shaping.block_reduce(image_data, downsample_factor, ridge=ridge)
    elif downsample_factor > 1:
        logger.info(f"Downsampling by factor {downsample_factor} (linear)...")
        has_blanks = not finite_mask.all()
        if has_blanks:
            # Interpolating across NaNs would spread each one into its
            # neighbours, so fill them first and carry the blank area along as
            # a coverage map, then blank again where coverage is below half.
            fill = nan_value if nan_value is not None else np.min(image_data[finite_mask])
            image_data[~finite_mask] = fill
            coverage = zoom(finite_mask.astype(np.float32), 1 / downsample_factor, order=1)
        image_data = zoom(image_data, 1 / downsample_factor, order=1)
        if image_data.size == 0:
            raise ValueError("Downsampling resulted in empty image.")
        if has_blanks:
            image_data[coverage < 0.5] = np.nan

    ny, nx = image_data.shape
    logger.info(f"Processing image dimensions: {nx}x{ny}")

    finite_mask = np.isfinite(image_data)
    if not np.all(finite_mask):
        num_non_finite = np.sum(~finite_mask)
        logger.warning(f"Found {num_non_finite} NaN/inf values.")
        if nan_value is not None:
            image_data[~finite_mask] = nan_value
        else:
            # Placeholder, will be replaced with min after clipping
            image_data[~finite_mask] = 0

    # Clipping
    if clip_percentile and 0 < clip_percentile < 50:
        logger.info(f"Clipping data to {clip_percentile}% - {100-clip_percentile}% percentile range.")
        data_for_percentile = image_data[finite_mask] if not np.all(finite_mask) else image_data
        if data_for_percentile.size > 0:
            min_val = np.percentile(data_for_percentile, clip_percentile)
            max_val = np.percentile(data_for_percentile, 100 - clip_percentile)
            image_data = np.clip(image_data, min_val, max_val)
            # Update NaNs if they weren't explicitly set
            if not np.all(finite_mask) and nan_value is None:
                image_data[~finite_mask] = min_val
        else:
            logger.warning("Cannot calculate percentiles, no finite data available.")
    elif not np.all(finite_mask) and nan_value is None:
        # If no clipping but has NaNs, use global min
        min_val = np.nanmin(image_data)
        image_data[~finite_mask] = min_val

    # Scaling
    #
    # log1p(x) and arcsinh(x) only compress once x is of order 1 or more, so on
    # data whose values are small in absolute terms - Jy/beam radio maps are
    # typically 1e-5 to 1e-1 - both are indistinguishable from a linear stretch.
    # `soften` divides through by a multiple of the noise before stretching,
    # which puts the data in units where the compression actually bites and
    # makes the result independent of the file's brightness units.
    if log_scale or asinh_scale:
        min_data = np.min(image_data)
        image_data = image_data - min_data
        if soften and soften > 0:
            scale = soften * noise
            if scale <= 0:
                logger.warning("Cannot soften: noise estimate is zero. Ignoring --soften.")
            else:
                logger.info(f"Softening the stretch by {soften:g} x rms ({scale:.4g}); "
                            f"peak is {np.max(image_data)/scale:.0f} softening units.")
                image_data = image_data / scale
        elif np.max(image_data) < 1.0:
            logger.warning(
                f"Data peaks at {np.max(image_data):.3g}, so log/arcsinh scaling will behave "
                "almost linearly and compress nothing. Use --soften N (in units of the noise "
                "rms) to get a real stretch."
            )

    if log_scale:
        logger.info("Applying log1p scaling...")
        image_data = np.log1p(image_data)
    elif asinh_scale:
        logger.info("Applying arcsinh scaling...")
        image_data = np.arcsinh(image_data)
    elif tactile:
        if tactile_edge is None:
            sky = float(np.median(image_data))
            tactile_edge = sky + 3 * robust_rms(image_data)
            logger.info(f"Tactile edge at sky + 3 x noise = {tactile_edge:.4g}")
        logger.info(f"Applying tactile stretch (edge {tactile_edge:.4g}, edge step {tactile_edge_step:g}, "
                    f"contrast {tactile_contrast:g})...")
        image_data = shaping.tactile_stretch(image_data, tactile_edge, peak_percentile=tactile_peak_percentile,
                                             edge_step=tactile_edge_step, contrast=tactile_contrast)

    if invert:
        logger.info("Inverting data height...")
        image_data = np.max(image_data) - image_data

    # Normalization to [0, 1] range
    min_val, max_val = np.min(image_data), np.max(image_data)
    data_range = max_val - min_val
    if data_range == 0:
        logger.warning("Data range is zero. Model will be flat.")
        normalized_data = np.zeros_like(image_data)
    else:
        normalized_data = (image_data - min_val) / data_range

    return normalized_data

def generate_mesh(
    normalized_data,
    longest_side_mm=None,
    max_height_mm=10.0,
    base_thickness_mm=2.0,
    smoothing_sigma=0,
    border_width_mm=0.0,
    border_height_mm=0.0,
    label_text="",
    label_depth_mm=None,
    label_font_size_mm=5.0,
    label_mode="text",
    label_area_mm=0.0,
    label_line_spacing_mm=None,
    image_labels=()
):
    """
    Generates a 3D mesh from the normalized image data.

    label_text is one line, several lines separated by newlines, or a list of
    lines. label_font_size_mm is the cap height of text labels.
    label_depth_mm is the height of the raised text, or of the braille dots in
    braille mode. If None, it defaults to DEFAULT_TEXT_DEPTH_MM for text and
    DEFAULT_BRAILLE_DOT_HEIGHT_MM for braille. label_line_spacing_mm is the
    distance between line centres (see generate_label_mesh).

    image_labels places extra label meshes on the image itself: a list of
    (mesh, x, y, z_mm), where the mesh is centred on x=0, y=0 and stands on
    z=0 (as from generate_label_mesh), (x, y) is its centre in pixels of
    normalized_data (row 0 at the bottom), and z_mm is the surface height
    it stands on.

    Raises ValueError if the label cannot be made, e.g. braille that does not
    fit across the model (braille cannot be scaled down and stay readable).
    """
    if label_mode not in ("text", "braille"):
        raise ValueError(f"label_mode must be 'text' or 'braille', not {label_mode!r}")

    ny, nx = normalized_data.shape

    # Calculate scaling factor (mm per pixel)
    if longest_side_mm and longest_side_mm > 0:
        scale_factor = longest_side_mm / max(nx, ny)
        logger.info(f"Scale factor: {scale_factor:.4f} mm/pixel")
    else:
        scale_factor = 1.0
        logger.info("Using 1 mm/pixel scaling (no longest side given)")
        if max(nx, ny) > 300:
            logger.warning(f"Without a longest side the model will be {max(nx, ny)} mm long; "
                           "set longest_side_mm (--longest-side) to fit your printer.")

    # Apply height scaling and base thickness
    z_data = base_thickness_mm + (normalized_data * max_height_mm)

    # Smoothing
    if smoothing_sigma > 0:
        logger.info(f"Applying Gaussian smoothing (sigma={smoothing_sigma} pixels)...")
        z_data = gaussian_filter(z_data, sigma=smoothing_sigma)

    bp = lp = 0   # border and label-strip pixels, to place image labels
    # Add border
    if border_width_mm > 0:
        border_pixels = int(round(border_width_mm / scale_factor))
        if border_pixels > 0:
            logger.info(f"Adding border: {border_width_mm}mm ({border_pixels} pixels)")
            ny_new = ny + 2 * border_pixels
            nx_new = nx + 2 * border_pixels
            # Border height is relative to the base thickness
            z_bordered = np.full((ny_new, nx_new), base_thickness_mm + border_height_mm, dtype=np.float32)
            z_bordered[border_pixels:border_pixels+ny, border_pixels:border_pixels+nx] = z_data
            z_data = z_bordered
            ny, nx = ny_new, nx_new
            bp = border_pixels
        else:
            logger.warning("Border width too small to represent at this scale, skipping.")

    # Add extra flat area for the label at the bottom
    if label_area_mm > 0:
        label_area_pixels = int(round(label_area_mm / scale_factor))
        if label_area_pixels > 0:
            logger.info(f"Adding dedicated label area: {label_area_mm}mm ({label_area_pixels} pixels)")
            # Prepend flat rows to the start of z_data (corresponds to negative Y)
            ny_new = ny + label_area_pixels
            z_label_area = np.full((ny_new, nx), base_thickness_mm + border_height_mm, dtype=np.float32)
            z_label_area[label_area_pixels:, :] = z_data
            z_data = z_label_area
            ny = ny_new
            lp = label_area_pixels

    # Create grid and center it around (0,0)
    x = (np.arange(nx, dtype=np.float32) - (nx - 1) / 2.0) * scale_factor
    y = (np.arange(ny, dtype=np.float32) - (ny - 1) / 2.0) * scale_factor
    xx, yy = np.meshgrid(x, y)

    num_vertices = nx * ny
    vertices = np.zeros((num_vertices * 2, 3), dtype=np.float32)
    
    # Top surface vertices
    vertices[:num_vertices, 0] = xx.flatten()
    vertices[:num_vertices, 1] = yy.flatten()
    vertices[:num_vertices, 2] = z_data.flatten()
    
    # Bottom surface vertices
    vertices[num_vertices:, 0] = xx.flatten()
    vertices[num_vertices:, 1] = yy.flatten()
    vertices[num_vertices:, 2] = 0.0

    logger.info("Generating mesh faces (vectorized)...")
    
    # Indices for grid cells
    i, j = np.meshgrid(np.arange(nx - 1), np.arange(ny - 1))
    v00 = (j * nx + i).flatten()
    v10 = (j * nx + (i + 1)).flatten()
    v01 = ((j + 1) * nx + i).flatten()
    v11 = ((j + 1) * nx + (i + 1)).flatten()

    # Top surface faces (pointing up)
    top_faces_1 = np.column_stack([v00, v10, v01])
    top_faces_2 = np.column_stack([v10, v11, v01])
    top_faces = np.vstack([top_faces_1, top_faces_2])

    # Bottom surface faces (pointing down)
    v00_b = v00 + num_vertices
    v10_b = v10 + num_vertices
    v01_b = v01 + num_vertices
    v11_b = v11 + num_vertices
    
    bottom_faces_1 = np.column_stack([v00_b, v01_b, v10_b])
    bottom_faces_2 = np.column_stack([v10_b, v01_b, v11_b])
    bottom_faces = np.vstack([bottom_faces_1, bottom_faces_2])

    # Side faces
    side_faces_list = []
    
    # Y-edges (top and bottom in image coordinates)
    for i_edge in [0, ny - 1]:
        idx = np.arange(nx - 1)
        v_t0 = i_edge * nx + idx
        v_t1 = i_edge * nx + idx + 1
        v_b0 = v_t0 + num_vertices
        v_b1 = v_t1 + num_vertices
        
        if i_edge == 0:
            side_faces_list.append(np.column_stack([v_t0, v_b0, v_b1]))
            side_faces_list.append(np.column_stack([v_t0, v_b1, v_t1]))
        else:
            side_faces_list.append(np.column_stack([v_t0, v_b1, v_b0]))
            side_faces_list.append(np.column_stack([v_t0, v_t1, v_b1]))

    # X-edges (left and right in image coordinates)
    for j_edge in [0, nx - 1]:
        idx = np.arange(ny - 1)
        v_t0 = idx * nx + j_edge
        v_t1 = (idx + 1) * nx + j_edge
        v_b0 = v_t0 + num_vertices
        v_b1 = v_t1 + num_vertices
        
        if j_edge == 0:
            side_faces_list.append(np.column_stack([v_t0, v_b1, v_b0]))
            side_faces_list.append(np.column_stack([v_t0, v_t1, v_b1]))
        else:
            side_faces_list.append(np.column_stack([v_t0, v_b0, v_b1]))
            side_faces_list.append(np.column_stack([v_t0, v_b1, v_t1]))

    all_faces = np.vstack([top_faces, bottom_faces] + side_faces_list)
    
    logger.info(f"Mesh generated: {len(vertices)} vertices, {len(all_faces)} faces")
    
    stl_mesh = mesh.Mesh(np.zeros(all_faces.shape[0], dtype=mesh.Mesh.dtype))
    stl_mesh.vectors = vertices[all_faces]
    
    # Add the label (one or more lines of raised text or braille) if requested.
    # Errors here propagate: silently saving a model without its label would go
    # unnoticed in a batch build.
    lines = _label_lines(label_text)
    if lines:
        logger.info(f"Adding {label_mode} label: {' / '.join(lines)!r}")
        model_width = (nx - 1) * scale_factor
        # Keep the label off the border, and at least 5% in from each edge.
        margin = max(border_width_mm, 0.05 * model_width)
        text_mesh, block_h = generate_label_mesh(
            lines, mode=label_mode, depth_mm=label_depth_mm, text_height_mm=label_font_size_mm,
            line_spacing_mm=label_line_spacing_mm, max_width_mm=model_width - 2 * margin)

        bottom_edge_y = -(ny - 1) / 2.0 * scale_factor
        if label_area_mm > 0:
            # Centre the block of lines in the dedicated label strip.
            if block_h > label_area_mm:
                logger.warning(f"Label is {block_h:.1f} mm tall but the label area is only "
                               f"{label_area_mm:g} mm; increase --label-area.")
            ty = bottom_edge_y + label_area_mm / 2.0
        else:
            # No strip: sit the block just above the bottom edge.
            ty = bottom_edge_y + 0.6 * block_h

        tz = base_thickness_mm + border_height_mm
        text_mesh.translate([0, ty, tz], inplace=True)

        # 3. Convert pyvista mesh to numpy-stl
        tri_text = text_mesh.triangulate()
        text_faces = tri_text.faces.reshape(-1, 4)[:, 1:]
        text_vertices = tri_text.points

        text_stl_data = np.zeros(text_faces.shape[0], dtype=mesh.Mesh.dtype)
        text_stl = mesh.Mesh(text_stl_data)
        text_stl.vectors = text_vertices[text_faces]

        # 4. Merge meshes
        stl_mesh = mesh.Mesh(np.concatenate([stl_mesh.data, text_stl.data]))
        logger.info(f"Label added successfully ({len(text_faces)} faces)")

    # Labels on the image itself, at image pixel positions.
    for lab_mesh, px, py, z_mm in image_labels or ():
        m = lab_mesh.copy()
        m.translate([(px + bp - (nx - 1) / 2.0) * scale_factor,
                     (py + bp + lp - (ny - 1) / 2.0) * scale_factor, z_mm], inplace=True)
        tri = m.triangulate()
        faces = tri.faces.reshape(-1, 4)[:, 1:]
        extra = mesh.Mesh(np.zeros(len(faces), dtype=mesh.Mesh.dtype))
        extra.vectors = tri.points[faces]
        stl_mesh = mesh.Mesh(np.concatenate([stl_mesh.data, extra.data]))

    return stl_mesh

def main():
    parser = argparse.ArgumentParser(
        description="AstroTouch: Convert astronomical FITS images to 3D printable STL files.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        allow_abbrev=False
    )
    
    # Files
    parser.add_argument("input", type=Path, help="Input FITS file, or an image (PNG, JPEG, TIFF)")
    tpl = parser.add_argument_group(
        "templates",
        "Start from a plaque format and/or an image recipe; any option you also give overrides it. "
        "Run with --list-templates to see them.")
    tpl.add_argument("--plaque", metavar="NAME", help="Plaque format: large, standard or small")
    tpl.add_argument("--recipe", metavar="NAME",
                     help="Image recipe: plain, radio, radio_filaments, galaxy_smooth, galaxy_knotty, nebula, "
                          "press_picture")
    tpl.add_argument("--template", type=Path, metavar="JSON",
                     help="A template saved from the GUI (or its .astrotouch.json)")
    tpl.add_argument("--list-templates", action="store_true", help="Describe the templates and exit")
    parser.add_argument("output", type=Path, nargs='?', help="Output STL file (default: input_name.stl)")
    
    # Selection
    parser.add_argument("--hdu", type=int, help="HDU index (default: auto-detect first 2D image)")
    
    # Physical dimensions
    parser.add_argument("--longest-side", type=float, default=150.0, metavar="MM",
                        help="Length of the image's longest side in mm (0 for 1 mm per pixel)")
    parser.add_argument("--max-height", type=float, default=10.0, metavar="MM", help="Maximum feature height above base")
    parser.add_argument("--base-thickness", type=float, default=2.0, metavar="MM", help="Thickness of the base plate")
    
    # Processing
    parser.add_argument("--clip", type=float, default=None, metavar="PERCENT",
                        help="Percentile to clip outliers (0 to disable; default 1, or 0 with --clean or --tactile)")
    parser.add_argument("--smooth", type=float, default=1.0, metavar="SIGMA",
                        help="Gaussian smoothing sigma in model pixels (0 to disable)")
    parser.add_argument("--downsample", type=int, default=1, metavar="FACTOR", help="Downsample factor (e.g., 2 halves dimensions)")
    parser.add_argument("--downsample-method", choices=["mean", "linear"], default="mean",
                        help="Block mean (averages noise away) or linear interpolation (the old behaviour)")
    parser.add_argument("--ridge", type=float, default=0.0, metavar="F",
                        help="Mix this fraction of the block maximum into block-mean downsampling, keeping "
                             "thin features (filaments) continuous; use on cleaned maps only")
    parser.add_argument("--min-pixel", type=float, default=0.3, metavar="MM",
                        help="Downsample automatically so model pixels are at least this size; finer than "
                             "a 0.4 mm nozzle only bloats the STL (0 to disable)")
    parser.add_argument("--nan-value", type=float, help="Value to replace NaNs with (default: minimum)")
    parser.add_argument("--label", type=str, action="append", default=[], metavar="TEXT",
                        help="Label line to add at the bottom edge; repeat for more lines "
                             "(e.g. --label 'Centaurus A' --label 'Radio')")
    parser.add_argument("--label-spacing", type=float, default=None, metavar="MM",
                        help=f"Distance between label line centres (default: {BRAILLE_LINE_SPACING_MM:g} for "
                             f"braille, {TEXT_LINE_SPACING:g} x --label-size for text)")
    parser.add_argument("--label-mode", type=str, choices=["text", "braille"], default="text", help="Mode for the label (standard text or braille)")
    parser.add_argument("--label-depth", type=float, default=None, metavar="MM",
                        help=f"Height of the raised text, or of the braille dots in braille mode "
                             f"(default: {DEFAULT_TEXT_DEPTH_MM} for text, {DEFAULT_BRAILLE_DOT_HEIGHT_MM} for braille)")
    parser.add_argument("--label-size", type=float, default=5.0, metavar="MM", help="Cap height of text labels in mm")
    parser.add_argument("--label-area", type=float, default=0.0, help="Width of dedicated flat area at the bottom for the label (mm)")
    
    # Scaling options
    scale_group = parser.add_mutually_exclusive_group()
    scale_group.add_argument("--log", action="store_true", help="Apply log(1+x) scaling to pixels")
    scale_group.add_argument("--asinh", action="store_true", help="Apply arcsinh scaling (often better for astronomy)")
    scale_group.add_argument("--tactile", action="store_true",
                             help="Tactile stretch: log brightness from an outer edge step to a capped peak, "
                                  "with a local-contrast boost (see 'tactile stretch' options)")
    parser.add_argument("--soften", type=float, default=0.0, metavar="N",
                        help="Divide by N x the noise rms before --log/--asinh, so the stretch actually "
                             "compresses whatever the data's units (e.g. Jy/beam)")
    tg = parser.add_argument_group("tactile stretch")
    tg.add_argument("--edge", type=float, metavar="VALUE",
                    help="Outer isophote in data units: fainter is flat base (default: sky + 3 x noise; "
                         "with --clean, --detect-sigma)")
    tg.add_argument("--edge-step", type=float, default=0.2, metavar="F",
                    help="Height of the outline step, as a fraction of full height")
    tg.add_argument("--contrast", type=float, default=1.5, metavar="F",
                    help="Local-contrast boost: lifts arms, bars, knots and shells (0 for a plain log stretch)")
    tg.add_argument("--peak-percentile", type=float, default=99.9, metavar="P",
                    help="Cap the brightness here so a bright core does not take the whole height range")
    
    parser.add_argument("--invert", action="store_true", help="Invert height map (brightest pixels become lowest)")
    
    # Border
    parser.add_argument("--border-width", type=float, default=0.0, metavar="MM", help="Border width in mm")
    parser.add_argument("--border-height", type=float, default=0.0, metavar="MM", help="Border height above base in mm")
    
    # Noise cleaning
    cg = parser.add_argument_group(
        "noise cleaning",
        "Mask the noise and flatten empty sky before the stretch (for radio maps and other "
        "noisy images). The cleaned map is in units of the noise rms, so combine with --asinh.")
    cg.add_argument("--clean", action="store_true", help="Enable noise cleaning")
    cg.add_argument("--beam", type=float, nargs=2, metavar=("BMAJ", "BMIN"),
                    help="Beam FWHM in pixels (default: from the BMAJ/BMIN header keywords)")
    cg.add_argument("--detect-sigma", type=float, default=3.0, metavar="N", help="Island boundary, in sigma")
    cg.add_argument("--peak-sigma", type=float, default=5.0, metavar="N",
                    help="An island must peak above this (at any scale) to be kept")
    cg.add_argument("--min-island", type=float, metavar="PX",
                    help="Smallest island kept, in pixels (default: half a beam area, or 12)")
    cg.add_argument("--smooth-scales", type=float, nargs="+", metavar="PX",
                    help="Gaussian sigmas of the extra detection scales (default: one beam sigma)")
    cg.add_argument("--noise-box", type=int, metavar="PX",
                    help="Measure a local noise map in boxes of this size (for mosaics with varying noise)")
    cg.add_argument("--save-cleaned", type=Path, metavar="FITS", help="Also write the cleaned map to this FITS file")

    # Source removal
    sg = parser.add_argument_group(
        "source removal",
        "Remove foreground stars and other compact sources (optical/IR images) and fill them by "
        "harmonic inpainting. Pixel positions are 0-based (x, y) in the input image; for PNG/JPEG, "
        "as shown in an image viewer (y down).")
    sg.add_argument("--mask-stars", action="store_true",
                    help="Mask foreground stars from Gaia DR3 (needs a celestial WCS and the network)")
    sg.add_argument("--gaia-cache", type=Path, metavar="CSV", help="Cache the Gaia star list in this file")
    sg.add_argument("--remove-points", action="store_true",
                    help="Remove compact sources found against a median filter (faint stars, clusters, "
                         "small background galaxies); skip for objects whose knots matter")
    sg.add_argument("--max-point-area", type=float, default=60, metavar="PX",
                    help="Largest compact source removed, in pixels")
    sg.add_argument("--point-sigma", type=float, default=5.0, metavar="N", help="Compact-source detection threshold")
    sg.add_argument("--point-box", type=int, default=9, metavar="PX",
                    help="Median filter width for compact sources: a few star widths (e.g. 21 for JWST pictures)")
    sg.add_argument("--mask", type=float, nargs=3, action="append", default=[], metavar=("X", "Y", "R"),
                    help="Also mask and inpaint this disc (pixels); repeatable")
    sg.add_argument("--protect", type=float, nargs=3, action="append", default=[], metavar=("X", "Y", "R"),
                    help="Never mask inside this disc (e.g. a galaxy nucleus); repeatable")

    # Compact islands and markers
    mg = parser.add_argument_group("compact islands and markers")
    mg.add_argument("--remove-compact", type=float, metavar="MM",
                    help="Remove isolated islands smaller than this across on the model (needs a flat "
                         "background: --clean or --tactile); they print as pips like braille dots")
    mg.add_argument("--compact-fraction", type=float, default=0.0, metavar="F",
                    help="Keep compact islands as low pips at this fraction of their height instead")
    mg.add_argument("--compact-cap", type=float, default=0.15, metavar="F",
                    help="...capped at this fraction of full height (with --compact-fraction)")
    mg.add_argument("--marker", type=float, nargs=2, action="append", default=[], metavar=("X", "Y"),
                    help="Add a marker dome at this pixel position; repeatable")
    mg.add_argument("--marker-radec", type=float, nargs=2, action="append", default=[], metavar=("RA", "DEC"),
                    help="Add a marker dome at this sky position in degrees (needs a WCS); repeatable")
    mg.add_argument("--marker-size", type=float, default=3.0, metavar="MM", help="Marker dome diameter")
    mg.add_argument("--marker-height", type=float, default=0.3, metavar="F",
                    help="Marker height above the surrounding relief, as a fraction of full height")

    # Hand-painted layers
    pg = parser.add_argument_group(
        "painted layers",
        "Masks painted in the GUI (saved next to the STL) or any FITS/PNG mask, stretched to cover the image.")
    pg.add_argument("--paint-out", type=Path, metavar="MASK", help="Paint out the non-zero areas of this mask")
    pg.add_argument("--paint-out-mode", choices=["fill", "flatten"], default="fill",
                    help="fill: inpaint from the surroundings (before cleaning); flatten: set to flat base")
    pg.add_argument("--protect-mask", type=Path, metavar="MASK",
                    help="Never remove stars or compact sources inside the non-zero areas of this mask")
    pg.add_argument("--texture-map", type=Path, metavar="MAP",
                    help="Raise a tactile texture where this map is 1 (dots), 2 (lines) or 3 (grid)")
    pg.add_argument("--texture-spacing", type=float, default=4.0, metavar="MM", help="Texture pattern spacing")
    pg.add_argument("--texture-height", type=float, default=0.6, metavar="MM", help="Texture height")

    # Shapes and labels on the image
    ag = parser.add_argument_group(
        "shapes and labels on the image",
        "Raised lines, arrows and rings, and braille or text labels placed on the image itself. "
        "Drawn in the GUI, or written as JSON (see the README).")
    ag.add_argument("--annotations", type=Path, metavar="JSON",
                    help="A JSON list of shapes and labels, or a GUI-saved .astrotouch.json")
    ag.add_argument("--shape-width", type=float, default=1.2, metavar="MM", help="Width of raised lines and rings")
    ag.add_argument("--shape-height", type=float, default=0.8, metavar="MM", help="Height of raised lines and rings")
    ag.add_argument("--image-label-mode", choices=["braille", "text"], default="braille",
                    help="Style of labels on the image")
    ag.add_argument("--image-label-pad", choices=["base", "raised"], default="base",
                    help="Labels on the image sit on a flat pad: flattened to the base, or a raised platform")
    ag.add_argument("--image-label-margin", type=float, default=2.0, metavar="MM",
                    help="Flat margin around each label on the image")

    # Misc
    parser.add_argument("--bed", type=float, nargs=2, metavar=("W", "D"),
                        help="Warn if the model does not fit a print bed of this size (mm)")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose debug output")

    if "--list-templates" in sys.argv[1:]:
        try:
            from . import pipeline as _p
        except ImportError:
            import pipeline as _p
        for heading, table in (("Plaque formats (--plaque)", _p.PLAQUE_FORMATS), ("Image recipes (--recipe)", _p.RECIPES)):
            print(heading)
            for name, t in table.items():
                print(f"  {name:16s} {t['title']}\n  {'':16s} {t['about']}\n")
        return
    args = parser.parse_args()

    # Options given explicitly on the command line (they override templates).
    explicit = set()
    for tok in sys.argv[1:]:
        if tok.startswith("--"):
            act = parser._option_string_actions.get(tok.split("=", 1)[0])
            if act is not None:
                explicit.add(act.dest)

    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
    if args.verbose:
        logger.setLevel(logging.DEBUG)

    # Determine output path if not provided
    if not args.output:
        args.output = args.input.with_suffix('.stl')
        logger.info(f"No output path provided, using: {args.output}")

    try:
        try:
            from . import pipeline
        except ImportError:
            import pipeline

        # 1. Load data
        builder = pipeline.ModelBuilder.from_file(args.input, args.hdu)
        header = builder.header
        is_picture = Path(args.input).suffix.lower() in IMAGE_SUFFIXES
        full_ny = builder.raw.shape[0]

        def to_array_xy(x, y):
            # Pictures are flipped on loading so row 0 is the bottom; positions
            # given as in an image viewer (y down) must be flipped to match.
            return (x, full_ny - 1 - y) if is_picture else (x, y)

        markers = [to_array_xy(x, y) for x, y in args.marker]
        if args.marker_radec:
            from astropy.wcs import WCS
            wcs = WCS(header).celestial if header is not None else None
            if wcs is None or not wcs.has_celestial:
                raise ValueError("--marker-radec needs a celestial WCS in the FITS header.")
            markers += [tuple(map(float, wcs.world_to_pixel_values(ra, dec))) for ra, dec in args.marker_radec]

        annotations = ()
        if args.annotations:
            import json
            loaded = json.loads(Path(args.annotations).read_text())
            if isinstance(loaded, dict):    # a GUI-saved .astrotouch.json
                loaded = loaded.get("settings", {}).get("annotations", [])
            annotations = tuple(loaded)
            logger.info(f"{len(annotations)} shape(s) and label(s) on the image")

        stretch = "log" if args.log else "asinh" if args.asinh else "tactile" if args.tactile else "linear"
        settings = pipeline.Settings(
            longest_side_mm=args.longest_side, max_height_mm=args.max_height,
            base_thickness_mm=args.base_thickness, clip=args.clip, smooth=args.smooth,
            downsample=args.downsample, downsample_method=args.downsample_method, ridge=args.ridge,
            min_pixel_mm=args.min_pixel, nan_value=args.nan_value, stretch=stretch, soften=args.soften,
            invert=args.invert, edge=args.edge, edge_step=args.edge_step, contrast=args.contrast,
            peak_percentile=args.peak_percentile, clean=args.clean,
            beam=tuple(args.beam) if args.beam else None, detect_sigma=args.detect_sigma,
            peak_sigma=args.peak_sigma, min_island=args.min_island,
            smooth_scales=tuple(args.smooth_scales) if args.smooth_scales else None, noise_box=args.noise_box,
            mask_stars=args.mask_stars, gaia_cache=str(args.gaia_cache) if args.gaia_cache else None,
            remove_points=args.remove_points, point_box=args.point_box, point_sigma=args.point_sigma,
            max_point_area=args.max_point_area,
            masks=tuple((*to_array_xy(x, y), r) for x, y, r in args.mask),
            protect=tuple((*to_array_xy(x, y), r) for x, y, r in args.protect),
            paint_mode=args.paint_out_mode, remove_compact=args.remove_compact,
            compact_fraction=args.compact_fraction, compact_cap=args.compact_cap, markers=tuple(markers),
            marker_size_mm=args.marker_size, marker_height=args.marker_height,
            texture_spacing_mm=args.texture_spacing, texture_height_mm=args.texture_height,
            border_width_mm=args.border_width, border_height_mm=args.border_height, label=tuple(args.label),
            label_mode=args.label_mode, label_depth_mm=args.label_depth, label_size_mm=args.label_size,
            label_area_mm=args.label_area, label_spacing_mm=args.label_spacing,
            annotations=annotations, shape_width_mm=args.shape_width, shape_height_mm=args.shape_height,
            image_label_mode=args.image_label_mode, image_label_pad=args.image_label_pad,
            image_label_margin_mm=args.image_label_margin)

        if args.plaque or args.recipe or args.template:
            custom = None
            if args.template:
                import json
                loaded = json.loads(Path(args.template).read_text())
                custom = loaded.get("settings", loaded) if isinstance(loaded, dict) else None
                if not isinstance(custom, dict):
                    raise ValueError(f"{args.template} is not a template file.")
                custom = {k: v for k, v in custom.items() if k not in pipeline.TEMPLATE_EXCLUDED}
            templ = pipeline.apply_template(pipeline.Settings(), plaque=args.plaque, recipe=args.recipe,
                                            custom=custom)
            # Settings fields set by options the user gave explicitly keep the user's values.
            dest_fields = {
                "longest_side": "longest_side_mm", "max_height": "max_height_mm",
                "base_thickness": "base_thickness_mm", "clip": "clip", "smooth": "smooth",
                "downsample": "downsample", "downsample_method": "downsample_method", "ridge": "ridge",
                "min_pixel": "min_pixel_mm", "nan_value": "nan_value", "log": "stretch", "asinh": "stretch",
                "tactile": "stretch", "soften": "soften", "invert": "invert", "edge": "edge",
                "edge_step": "edge_step", "contrast": "contrast", "peak_percentile": "peak_percentile",
                "clean": "clean", "beam": "beam", "detect_sigma": "detect_sigma", "peak_sigma": "peak_sigma",
                "min_island": "min_island", "smooth_scales": "smooth_scales", "noise_box": "noise_box",
                "mask_stars": "mask_stars", "remove_points": "remove_points", "point_box": "point_box",
                "point_sigma": "point_sigma", "max_point_area": "max_point_area",
                "remove_compact": "remove_compact", "compact_fraction": "compact_fraction",
                "compact_cap": "compact_cap", "marker_size": "marker_size_mm", "marker_height": "marker_height",
                "texture_spacing": "texture_spacing_mm", "texture_height": "texture_height_mm",
                "paint_out_mode": "paint_mode", "shape_width": "shape_width_mm",
                "shape_height": "shape_height_mm", "image_label_mode": "image_label_mode",
                "image_label_pad": "image_label_pad", "image_label_margin": "image_label_margin_mm",
                "border_width": "border_width_mm", "border_height": "border_height_mm",
                "label_mode": "label_mode", "label_depth": "label_depth_mm", "label_size": "label_size_mm",
                "label_area": "label_area_mm", "label_spacing": "label_spacing_mm"}
            mine = {dest_fields[d] for d in explicit if d in dest_fields}
            merged = templ.to_dict()
            merged.update({k: v for k, v in settings.to_dict().items()
                           if k in mine or k in pipeline.TEMPLATE_EXCLUDED})
            settings = pipeline.Settings(**{k: (tuple(v) if isinstance(v, list) else v) for k, v in merged.items()})
            used = ", ".join(x for x in (args.plaque and f"plaque '{args.plaque}'",
                                         args.recipe and f"recipe '{args.recipe}'",
                                         args.template and f"template {args.template.name}") if x)
            logger.info(f"Using {used}" + (f"; overridden by {', '.join(sorted(mine))}" if mine else ""))

        # Hand-painted layers (from the GUI, or any mask image)
        layers = pipeline.Layers(
            paint=pipeline.load_layer(args.paint_out, builder.raw.shape) if args.paint_out else None,
            protect=pipeline.load_layer(args.protect_mask, builder.raw.shape) if args.protect_mask else None,
            texture=pipeline.load_layer(args.texture_map, builder.raw.shape, integer=True)
            if args.texture_map else None)

        # 2. Height map: sources, cleaning, stretch, compact islands, markers, texture
        result = builder.height_map(settings, layers)
        if result.factor > args.downsample:
            logger.info(f"Downsampled by {result.factor} so model pixels are at least {args.min_pixel:g} mm "
                        "(--min-pixel 0 to disable).")
        if args.clean and args.save_cleaned:
            hdr = header.copy() if header is not None else fits.Header()
            hdr["BUNIT"] = "sigma"
            for line in result.clean.history:
                hdr["HISTORY"] = line
            fits.writeto(args.save_cleaned, result.clean.image, hdr, overwrite=True)
            logger.info(f"Cleaned map written to {args.save_cleaned}")
        if args.clean and stretch in ("linear", "log") and not args.soften:
            logger.warning("Cleaned maps are in units of the noise; --asinh or --tactile usually gives "
                           "faint emission much more height than a linear stretch.")
        if markers:
            logger.info(f"Added {len(markers)} marker dome(s); markers are not data, say so when presenting.")

        # 3. Generate the 3D mesh
        stl_mesh = builder.mesh(settings, result)
        pts = stl_mesh.vectors.reshape(-1, 3)
        size = (float(np.ptp(pts[:, 0])), float(np.ptp(pts[:, 1])), float(pts[:, 2].max()))
        logger.info(f"Model size: {size[0]:.0f} x {size[1]:.0f} x {size[2]:.1f} mm")
        if args.bed and not pipeline.fits_bed(size, args.bed):
            logger.warning(f"The model does not fit a {args.bed[0]:g} x {args.bed[1]:g} mm bed, "
                           "even rotated; reduce --longest-side.")

        # 4. Save the STL file
        logger.info(f"Saving STL to: {args.output}")
        stl_mesh.save(args.output)
        logger.info("Success! 3D model generated.")

    except KeyboardInterrupt:
        logger.info("\nInterrupted by user.")
        sys.exit(0)
    except Exception as e:
        logger.error(f"Error: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
