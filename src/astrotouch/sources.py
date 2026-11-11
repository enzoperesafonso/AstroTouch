"""
Removing stars and other compact sources from optical and infrared images.

Printed straight, foreground stars become sharp spikes that can't be told apart
by touch from star-forming knots; bright stars bring bleed trails and
diffraction spikes that print as straight ridges; faint stars, globular
clusters and background galaxies cover smooth galaxies in bumps.

- Foreground stars come from Gaia DR3 (network; the list can be cached). Only
  sources with significant parallax or proper motion and a clean single-star
  fit (RUWE < 1.4), or brighter than G = 15.5, count as stars, so clusters and
  knots inside a galaxy (which Gaia also lists) are left alone.
- Other compact sources are found as peaks well above a median filter.
- Masked pixels are filled by harmonic inpainting, which cannot create bumps or
  dips of its own.

Ported from bvi-models galaxy-morphology/clean_galaxy.py.
"""
import logging
import os
import urllib.error
import urllib.parse
import urllib.request

import numpy as np
from scipy import ndimage

logger = logging.getLogger(__name__)

# Gaia DR3 TAP services, tried in order: ESA's archive, then the ARI Heidelberg mirror.
GAIA_SERVICES = [
    ("ESA Gaia archive", "https://gea.esac.esa.int/tap-server/tap/sync"),
    ("ARI Heidelberg mirror", "https://gaia.ari.uni-heidelberg.de/tap/sync"),
]
GAIA_COLUMNS = ("ra", "dec", "phot_g_mean_mag", "parallax_over_error", "pmra", "pmdec",
                "pmra_error", "pmdec_error", "ruwe")
_FAILED = {}               # query -> (time, message): don't hammer dead services from a live preview
FAILURE_MEMORY_S = 60
BRIGHT_G = 15.5


def _mad(x):
    return 1.4826 * np.median(np.abs(x - np.median(x)))


def _query_gaia(query, timeout=120, retries=1, wait_s=3.0):
    """Runs a synchronous ADQL query on each Gaia service in turn; returns CSV bytes."""
    import time
    body = urllib.parse.urlencode(dict(REQUEST="doQuery", LANG="ADQL", FORMAT="csv", QUERY=query)).encode()
    problems = []
    for name, url in GAIA_SERVICES:
        for attempt in range(retries + 1):
            try:
                logger.info(f"Gaia: querying the {name}...")
                req = urllib.request.Request(url, body, headers={"User-Agent": "astrotouch"})
                with urllib.request.urlopen(req, timeout=timeout) as r:
                    raw = r.read()
                header = raw.split(b"\n", 1)[0].decode("utf-8", "replace").strip().split(",")
                if [c.strip().strip('"') for c in header] != list(GAIA_COLUMNS):
                    raise ValueError("unexpected response (not a Gaia table)")
                return raw
            except urllib.error.HTTPError as e:
                detail = e.read().decode("utf-8", "replace")
                msg = next((line.strip() for line in detail.splitlines() if "Message:" in line), "")
                problems.append(f"{name}: HTTP {e.code}" + (f" ({msg.split('Message:')[-1].strip()})" if msg else ""))
                transient = e.code in (429, 500, 502, 503, 504)
            except (urllib.error.URLError, TimeoutError, OSError) as e:
                problems.append(f"{name}: {getattr(e, 'reason', e)}")
                transient = True
            except ValueError as e:
                problems.append(f"{name}: {e}")
                transient = False
            if not transient or attempt == retries:
                break
            time.sleep(wait_s)
    raise RuntimeError("Could not get the Gaia star list. " + "; ".join(problems) +
                       ". Try again later, or pass a list you already have with --gaia-cache "
                       "(columns: " + ", ".join(GAIA_COLUMNS) + ").")


def gaia_stars(ra, dec, radius_deg, cache=None, mag_limit=20.5):
    """
    Foreground stars from Gaia DR3 within radius_deg of (ra, dec), as a
    structured array with ra, dec and phot_g_mean_mag. If `cache` is a path,
    the query result is stored there and reused. The ESA archive is tried
    first, then the ARI Heidelberg mirror.
    """
    import time
    if cache and os.path.exists(cache):
        logger.info(f"Gaia: using cached list {cache}")
        with open(cache, "rb") as f:
            raw = f.read()
    else:
        q = ("SELECT " + ", ".join(GAIA_COLUMNS) + " FROM gaiadr3.gaia_source WHERE 1=CONTAINS("
             f"POINT('ICRS', ra, dec), CIRCLE('ICRS', {ra:.6f}, {dec:.6f}, {radius_deg:.6f})) "
             f"AND phot_g_mean_mag < {mag_limit}")
        failed = _FAILED.get(q)
        if failed and time.time() - failed[0] < FAILURE_MEMORY_S:
            raise RuntimeError(failed[1])
        try:
            raw = _query_gaia(q)
        except RuntimeError as e:
            _FAILED[q] = (time.time(), str(e))
            raise
        _FAILED.pop(q, None)
        if cache:
            with open(cache, "wb") as f:
                f.write(raw)
    lines = raw.decode().strip().splitlines()
    if len(lines) < 2:
        return np.zeros(0, dtype=[("ra", float), ("dec", float), ("phot_g_mean_mag", float)])
    from io import StringIO
    t = np.atleast_1d(np.genfromtxt(StringIO("\n".join(lines)), delimiter=",", names=True))
    plx = np.nan_to_num(t["parallax_over_error"])
    pm = np.hypot(np.nan_to_num(t["pmra"] / t["pmra_error"]),
                  np.nan_to_num(t["pmdec"] / t["pmdec_error"]))
    clean_fit = np.nan_to_num(t["ruwe"], nan=9) < 1.4
    star = (((plx > 3) | (pm > 5)) & clean_fit) | (t["phot_g_mean_mag"] < BRIGHT_G)
    logger.info(f"Gaia: {int(star.sum())} foreground stars, {int((~star).sum())} other sources left alone.")
    return t[star]


def _protected(shape, protect, protect_mask=None):
    yy, xx = np.indices(shape)
    keep = np.zeros(shape, bool) if protect_mask is None else np.asarray(protect_mask, bool).copy()
    for x, y, r in protect or ():
        keep |= np.hypot(xx - x, yy - y) < r
    return keep


def star_mask(shape, wcs, stars, pixscale_arcsec, protect=(), protect_mask=None, counts=None):
    """
    Mask of the given stars: a disc growing with brightness (3" + 25" x
    10^(-0.2 (G - 10)); 1.8x larger below G = 9 for the scattered-light ring),
    plus a bleed trail and diffraction spikes for stars brighter than G = 12.
    `protect` is a list of (x, y, radius) pixel discs never masked (e.g. a
    galaxy nucleus or a nebula's central star).
    """
    ny, nx = shape
    yy, xx = np.mgrid[0:ny, 0:nx]
    mask = np.zeros(shape, bool)
    if len(stars) == 0:
        return mask
    x, y = wcs.world_to_pixel_values(stars["ra"], stars["dec"])
    for xi, yi, g in zip(x, y, stars["phot_g_mean_mag"]):
        # A "star" centred in a protected disc is the protected object itself
        # (Gaia lists bright galaxy nuclei): skip it entirely, not just the disc.
        if any(np.hypot(xi - px, yi - py) < pr for px, py, pr in protect or ()):
            continue
        if protect_mask is not None and 0 <= int(round(yi)) < ny and 0 <= int(round(xi)) < nx \
                and protect_mask[int(round(yi)), int(round(xi))]:
            continue
        r = np.clip(3 + 25 * 10 ** (-0.2 * (g - 10)), 3, 300) / pixscale_arcsec
        if g < 9:
            r *= 1.8
        r = max(r, 2.5)
        if not (-r < xi < nx + r and -r < yi < ny + r):
            continue
        if counts is not None:
            counts["stars"] = counts.get("stars", 0) + 1
        mask |= (xx - xi) ** 2 + (yy - yi) ** 2 < r * r
        if g < 12:
            L, w = 5 * r, max(2.5, 0.1 * r)
            mask |= (np.abs(yy - yi) < w) & (np.abs(xx - xi) < L)
            mask |= (np.abs(xx - xi) < w) & (np.abs(yy - yi) < L)
    return mask & ~_protected(shape, protect, protect_mask)


def point_sources(img, rms, max_area=60, box=9, nsig=5.0, protect=(), protect_mask=None):
    """
    Compact sources Gaia does not catch (faint stars, globular clusters, small
    background galaxies): islands more than nsig x the noise above a median
    filter `box` px wide, up to max_area px in area. Dilated by 2 px.
    """
    res = img - ndimage.median_filter(img, box)
    local = 1.4826 * ndimage.median_filter(np.abs(res), 31)
    det = res > nsig * np.maximum(local, rms)
    lab, n = ndimage.label(det)
    if n == 0:
        return np.zeros(img.shape, bool)
    area = ndimage.sum(det, lab, np.arange(1, n + 1))
    m = np.isin(lab, np.nonzero(area <= max_area)[0] + 1)
    m &= ~_protected(img.shape, protect, protect_mask)
    return ndimage.binary_dilation(m, iterations=2)


def inpaint(img, mask, iters=2000):
    """
    Harmonic inpainting: masked pixels solve Laplace's equation with the
    surrounding pixels as boundary values. A harmonic fill has no maxima or
    minima of its own, so a filled hole is never felt as a bump or a dip.
    """
    out = np.array(img, dtype=np.float64, copy=True)
    if not mask.any():
        return out
    good = (~mask).astype(float)
    todo = mask.copy()
    for s in (2, 4, 8, 16, 32, 64):
        w = ndimage.gaussian_filter(good, s)
        v = ndimage.gaussian_filter(np.where(mask, 0.0, out), s) / np.maximum(w, 1e-9)
        fill = todo & (w > 0.2)
        out[fill] = v[fill]
        todo &= ~fill
    out[todo] = np.median(out[~mask]) if (~mask).any() else 0.0
    ys, xs = np.nonzero(ndimage.binary_dilation(mask, iterations=1))
    sl = (slice(ys.min(), ys.max() + 1), slice(xs.min(), xs.max() + 1))
    sub, m = out[sl], mask[sl]
    for _ in range(iters):
        p = np.pad(sub, 1, mode="edge")
        avg = 0.25 * (p[:-2, 1:-1] + p[2:, 1:-1] + p[1:-1, :-2] + p[1:-1, 2:])
        sub[m] = avg[m]
    out[sl] = sub
    return out


def remove_sources(data, header=None, gaia=False, gaia_cache=None, points=False, point_area=60,
                   point_sigma=5.0, point_box=9, extra_masks=(), protect=(), extra_mask=None,
                   protect_mask=None, parts=False):
    """
    Removes foreground stars (Gaia; needs a celestial WCS in the header),
    compact sources (sharper than a median filter point_box px wide; widen it
    for images with broad stars), and hand-placed discs (extra_masks: (x, y, radius) in
    pixels) or painted (extra_mask, a boolean array), then fills them by
    harmonic inpainting. Nothing inside `protect` discs or `protect_mask` is
    masked, except what is painted out by hand. NaNs are kept as NaN.

    Returns (image, mask), or with parts=True (image, mask, parts): parts
    holds what was removed and why, for review: boolean masks "stars" (Gaia),
    "points" (compact sources), "discs" (extra_masks) and "painted", and the
    counts "n_stars" and "n_points".
    """
    img = np.array(data, dtype=np.float64, copy=True)
    valid = np.isfinite(img)
    fill = np.median(img[valid])
    img[~valid] = fill
    yy, xx = np.indices(img.shape)
    discs = np.zeros(img.shape, bool)
    for x, y, r in extra_masks or ():
        discs |= np.hypot(xx - x, yy - y) < r
    painted = np.zeros(img.shape, bool) if extra_mask is None else np.asarray(extra_mask, bool)
    star_m = np.zeros(img.shape, bool)
    point_m = np.zeros(img.shape, bool)
    counts = {"stars": 0}

    if gaia:
        from astropy.wcs import WCS
        from astropy.wcs.utils import proj_plane_pixel_scales
        if header is None:
            raise ValueError("Gaia star masking needs a FITS header with a celestial WCS.")
        wcs = WCS(header).celestial
        if not wcs.has_celestial:
            raise ValueError("Gaia star masking needs a celestial WCS in the FITS header.")
        ny, nx = img.shape
        pixscale = float(np.mean(proj_plane_pixel_scales(wcs))) * 3600
        ra, dec = wcs.pixel_to_world_values((nx - 1) / 2, (ny - 1) / 2)
        cra, cdec = wcs.pixel_to_world_values([0, nx - 1, 0, nx - 1], [0, 0, ny - 1, ny - 1])
        from astropy.coordinates import SkyCoord
        import astropy.units as u
        c0 = SkyCoord(float(ra) * u.deg, float(dec) * u.deg)
        radius = float(np.max(c0.separation(SkyCoord(cra * u.deg, cdec * u.deg)).deg)) * 1.05
        stars = gaia_stars(float(ra), float(dec), radius, cache=gaia_cache)
        star_m = star_mask(img.shape, wcs, stars, pixscale, protect=protect, protect_mask=protect_mask,
                           counts=counts)

    if points:
        border = max(2, int(0.08 * min(img.shape)))
        edge_px = np.concatenate([img[:border].ravel(), img[-border:].ravel(),
                                  img[:, :border].ravel(), img[:, -border:].ravel()])
        rms = _mad(edge_px)
        already = discs | star_m | painted
        point_m = point_sources(np.where(already, np.median(img), img), rms, max_area=point_area,
                                box=point_box, nsig=point_sigma, protect=protect, protect_mask=protect_mask)
        point_m &= ~already

    mask = (discs | painted | star_m | point_m) & valid
    n_points = int(ndimage.label(point_m)[1]) if points else 0
    found = []
    if gaia:
        found.append(f"{counts['stars']} Gaia star(s)")
    if points:
        found.append(f"{n_points} compact source(s)")
    logger.info(f"Source removal: {', '.join(found) + '; ' if found else ''}"
                f"inpainting {100 * mask.mean():.2f}% of the image.")
    out = inpaint(img, mask)
    out[~valid] = np.nan
    if parts:
        info = dict(stars=star_m & valid, points=point_m & valid, discs=discs & valid,
                    painted=painted & valid, n_stars=counts["stars"], n_points=n_points)
        return out.astype(np.float32), mask, info
    return out.astype(np.float32), mask
