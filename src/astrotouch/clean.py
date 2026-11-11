"""
Noise cleaning for tactile models of noisy images (radio maps in particular).

Printed straight from the FITS file, a noisy image makes a poor model:

- the noise prints as a rough, sandpaper texture over the whole plaque;
- noise peaks print as random bumps that can't be told apart from real sources;
- faint diffuse emission (radio lobes, supernova shells), often only 3-5 sigma
  per pixel, is lost in that texture;
- negative "bowls" around bright sources (interferometer artefacts) and blank
  survey edges print as dips.

clean_image() fixes this before the height stretch:

1. Noise: a robust rms (MAD of the quieter 90% of pixels), or optionally a
   local noise map measured in boxes, for mosaics whose noise varies.
2. Detection at several scales: the native map plus one or more
   Gaussian-smoothed copies (by default about one beam), in which diffuse
   emission rises well above its own noise.
3. Island mask: a pixel is emission if it is above `detect_sigma` at any scale;
   an island is kept only if it peaks above `peak_sigma` at some scale and is
   at least `min_area_px` in area.
4. Blend by brightness: faint emission comes from the first smoothed map, so it
   prints as a smooth surface; bright emission keeps native resolution.
5. Flat sky: everything outside the (slightly dilated, soft-edged) mask is
   exactly zero, so empty sky prints as a perfectly flat base and each source
   meets it on a short slope rather than a sheer wall.
6. Units of the noise rms: an asinh stretch (--asinh) then compresses from
   about 1 sigma upward, which gives faint emission real height. Use --clip 0
   afterwards so the true peaks are kept.

This generalises the cleaning used for the bvi-models MeerKAT sets
(Delhaize giant radio galaxies, MeerKAT Galactic Centre).
"""
from dataclasses import dataclass, field
import logging

import numpy as np
from scipy import ndimage

logger = logging.getLogger(__name__)

FWHM_PER_SIGMA = 2.0 * np.sqrt(2.0 * np.log(2.0))


def mad_std(values):
    """Standard deviation estimated from the median absolute deviation."""
    return 1.4826 * np.median(np.abs(values - np.median(values)))


def robust_rms(data, valid=None):
    """Noise rms from the MAD of the quieter 90% of pixels (ignores bright sources)."""
    vals = data[valid] if valid is not None else data[np.isfinite(data)]
    if vals.size == 0:
        raise ValueError("No valid pixels to estimate the noise from.")
    a = np.abs(vals)
    quiet = vals[a < np.percentile(a, 90)]
    if quiet.size == 0:
        quiet = vals
    return float(mad_std(quiet))


def beam_from_header(header):
    """
    The restoring beam (BMAJ, BMIN) in pixels from a FITS header, or None if
    the header lacks the beam or the pixel scale.
    """
    if header is None or "BMAJ" not in header or "BMIN" not in header:
        return None
    if "CDELT2" in header and header["CDELT2"]:
        pix = abs(float(header["CDELT2"]))
    elif all(k in header for k in ("CD1_1", "CD1_2", "CD2_1", "CD2_2")):
        pix = np.sqrt(abs(header["CD1_1"] * header["CD2_2"] - header["CD1_2"] * header["CD2_1"]))
    else:
        return None
    if pix == 0:
        return None
    return float(header["BMAJ"]) / pix, float(header["BMIN"]) / pix


# --- Local noise map (for mosaics) ----------------------------------------------

def _boxed(stat, a, valid, box):
    """Per-box statistic over valid pixels; NaN where a box is mostly invalid."""
    ny, nx = a.shape
    by, bx = max(1, ny // box), max(1, nx // box)
    grid = np.full((by, bx), np.nan)
    for j in range(by):
        for i in range(bx):
            sl = (slice(j * ny // by, (j + 1) * ny // by), slice(i * nx // bx, (i + 1) * nx // bx))
            t = a[sl][valid[sl]]
            if t.size > 0.3 * (a[sl].size):
                grid[j, i] = stat(t)
    return grid


def _grid_to_map(grid, shape, box):
    """Interpolate a box grid back to a smooth full-size map."""
    g = np.where(np.isfinite(grid), grid, np.nanmedian(grid))
    g = ndimage.median_filter(g, 3, mode="nearest")
    m = ndimage.zoom(g, (shape[0] / g.shape[0], shape[1] / g.shape[1]), order=1, mode="nearest")
    m = m[:shape[0], :shape[1]]
    if m.shape != shape:  # zoom can round one pixel short
        m = np.pad(m, [(0, shape[0] - m.shape[0]), (0, shape[1] - m.shape[1])], mode="edge")
    return ndimage.gaussian_filter(m, box / 2)


def local_noise(img, valid, box):
    """
    Local small-scale noise map, and the reference (quiet-sky) rms.

    The noise is measured on a high-pass residual (image minus a 7-pixel median
    filter), which ignores smooth emission and negative bowls, in boxes of
    `box` pixels. It is calibrated against the image MAD in the quietest 25%
    of boxes, and floored at the 10th percentile.

    Returns (rms_map, rms_ref, quiet_boxes, image_box_mads).
    """
    hp = img - ndimage.median_filter(img, 7)
    hp_grid = _boxed(mad_std, hp, valid, box)
    img_grid = _boxed(mad_std, img, valid, box)
    if not np.isfinite(hp_grid).any():
        raise ValueError(f"Noise box of {box} px is too large for a {img.shape[1]}x{img.shape[0]} image "
                         "with this much blank area.")
    quiet = hp_grid <= np.nanpercentile(hp_grid, 25)
    cal = np.nanmedian(img_grid[quiet] / hp_grid[quiet])
    floor = np.nanpercentile(hp_grid, 10) * cal
    rms_map = np.maximum(_grid_to_map(hp_grid * cal, img.shape, box), floor)
    rms_ref = float(np.nanmedian(img_grid[quiet]))
    return rms_map, rms_ref, quiet, img_grid


# --- Cleaning -----------------------------------------------------------------------

@dataclass
class CleanResult:
    image: np.ndarray            # cleaned map in units of rms_ref, zero background
    mask: np.ndarray             # emission kept (before dilation)
    rms: float                   # reference noise rms, in the input's units
    islands: int                 # islands detected
    kept: int                    # islands kept
    scales_px: tuple             # Gaussian sigmas of the extra detection scales
    min_area_px: float
    history: list = field(default_factory=list)  # lines for FITS HISTORY cards


def clean_image(
    data,
    beam_px=None,
    smooth_scales_px=None,
    detect_sigma=3.0,
    peak_sigma=5.0,
    min_area_px=None,
    blend_sigma=(4.0, 12.0),
    noise_box_px=None,
    dilate_px=2,
    edge_softening_px=1.5,
):
    """
    Masks noise and flattens empty sky; returns a CleanResult whose image is in
    units of the noise rms with an exactly zero background.

    beam_px: the beam FWHM in pixels, as (bmaj, bmin) or one number. It sets
        the default detection scales and minimum island area. Without it, 2 px
        is assumed.
    smooth_scales_px: Gaussian sigmas (pixels) of the extra detection scales.
        The first one supplies the smooth surface for faint emission. Default:
        one beam sigma (at least 1 px), as used for the Delhaize giant radio
        galaxies. Add broader scales to lift very faint, very extended emission
        (the MeerKAT Galactic Centre mosaic used 3 and 9 px with a ~1.8 px
        beam), at the cost of wider, softer island edges.
    detect_sigma, peak_sigma: island boundary, and the level an island must
        reach somewhere (at any scale) to be kept.
    min_area_px: smallest island kept. Default: half a beam area (or 12 px).
    blend_sigma: (lo, hi) S/N over which the map goes from smoothed to native.
    noise_box_px: if set, measure a local noise map in boxes of this size (for
        mosaics with varying noise) instead of one global rms.
    dilate_px, edge_softening_px: grow the mask slightly and taper it, so
        sources meet the base on a slope.
    """
    img = np.array(data, dtype=np.float64, copy=True)
    if img.ndim != 2:
        raise ValueError("clean_image needs a 2D image.")
    valid = np.isfinite(img)
    if not valid.any():
        raise ValueError("Image contains no finite pixels.")
    img[~valid] = 0.0

    if beam_px is None:
        bmaj = bmin = 2.0
        beam_known = False
    else:
        bmaj, bmin = (beam_px, beam_px) if np.isscalar(beam_px) else beam_px
        beam_known = True
    beam_sigma = np.sqrt(bmaj * bmin) / FWHM_PER_SIGMA
    beam_area = np.pi * bmaj * bmin / (4 * np.log(2))

    if smooth_scales_px is None:
        smooth_scales_px = (max(beam_sigma, 1.0),)
    smooth_scales_px = tuple(float(s) for s in np.atleast_1d(smooth_scales_px))
    if not smooth_scales_px or min(smooth_scales_px) <= 0:
        raise ValueError("smooth_scales_px must be one or more positive numbers.")
    if min_area_px is None:
        min_area_px = 0.5 * beam_area if beam_known else 12.0

    # 1. Noise, and S/N at each scale.
    smooths = [ndimage.gaussian_filter(img, s) for s in smooth_scales_px]
    if noise_box_px:
        rms_map, rms_ref, quiet, img_grid = local_noise(img, valid, int(noise_box_px))
        snrs = [img / rms_map]
        for sm in smooths:
            sm_grid = _boxed(mad_std, sm, valid, int(noise_box_px))
            f = np.nanmedian(sm_grid[quiet] / img_grid[quiet])   # noise reduction by smoothing
            snrs.append(sm / (rms_map * f))
        noise_desc = f"local noise in {int(noise_box_px)} px boxes, quiet-sky rms {rms_ref:.4g}"
    else:
        rms_ref = robust_rms(img, valid)
        if rms_ref > 0:
            snrs = [img / rms_ref] + [sm / robust_rms(sm, valid) for sm in smooths]
        noise_desc = f"global rms {rms_ref:.4g}"
    if not rms_ref > 0:
        raise ValueError("The noise estimate is zero: the image has no noise to clean "
                         "(is it already cleaned, or a display image?).")

    # 2-3. Multi-scale detection and island mask.
    detect = np.zeros(img.shape, bool)
    for snr in snrs:
        detect |= snr > detect_sigma
    detect &= valid
    labels, n = ndimage.label(detect)
    idx = np.arange(1, n + 1)
    if n:
        peak = np.max([ndimage.maximum(snr, labels, idx) for snr in snrs], axis=0)
        area = ndimage.sum(detect, labels, idx)
        keep = (peak > peak_sigma) & (area >= min_area_px)
    else:
        keep = np.zeros(0, bool)
    mask = np.isin(labels, idx[keep])
    if not mask.any():
        logger.warning("Cleaning kept no emission: the model will be flat. "
                       "Try a lower --peak-sigma or --detect-sigma.")

    # 4. Blend: faint emission from the first smoothed map, bright from native.
    lo, hi = blend_sigma
    w = np.clip((snrs[0] - lo) / (hi - lo), 0, 1)
    w = ndimage.gaussian_filter(w, 2)
    blended = (w * img + (1 - w) * smooths[0]) / rms_ref

    # 5. Flat sky with soft-edged sources.
    mask_d = ndimage.binary_dilation(mask, iterations=int(dilate_px)) if dilate_px > 0 else mask
    if edge_softening_px > 0:
        taper = np.clip(ndimage.gaussian_filter(mask_d.astype(float), edge_softening_px) * 1.5, 0, 1)
    else:
        taper = 1.0
    out = np.clip(blended, 0, None) * taper
    out[~mask_d | ~valid] = 0.0

    scales_txt = ", ".join(f"{s:g}" for s in smooth_scales_px)
    logger.info(f"Cleaning: {noise_desc}; detection scales native + sigma {scales_txt} px; "
                f"kept {int(keep.sum())} of {n} islands (>{detect_sigma:g} sigma, peak >{peak_sigma:g} "
                f"sigma, >= {min_area_px:.0f} px); {100 * mask_d.mean():.1f}% of the image is emission; "
                f"peak {out.max():.0f} sigma.")
    history = [
        "AstroTouch clean: multi-scale island mask, background zeroed, units of rms",
        f"  {noise_desc}",
        f"  scales native + {scales_txt} px; detect {detect_sigma:g}, peak {peak_sigma:g} sigma; "
        f"min area {min_area_px:.0f} px",
        f"  kept {int(keep.sum())} of {n} islands",
    ]
    return CleanResult(image=out.astype(np.float32), mask=mask, rms=rms_ref, islands=int(n),
                       kept=int(keep.sum()), scales_px=smooth_scales_px, min_area_px=float(min_area_px),
                       history=history)
