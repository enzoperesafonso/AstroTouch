"""Regression tests for the core engine: braille encoding, labels, preprocessing, mesh."""
import sys
from pathlib import Path

import numpy as np
import pytest
import pyvista as pv

from astrotouch import core  # noqa: E402
from astrotouch.core import ueb_cells  # noqa: E402

L = core._LETTERS
CAP, NUM, G1, SP = core._CAPITAL, core._NUMBER, core._GRADE1, core._SPACE


# --- UEB encoding -----------------------------------------------------------

@pytest.mark.parametrize("text, expected", [
    # Mixed-case word: a capital sign per capital, not the capitals-word
    # indicator (which would read "SBB").
    ("SBb", [CAP, L["s"], CAP, L["b"], L["b"]]),
    ("Sa", [CAP, L["s"], L["a"]]),
    # All-capitals word gets the capitals-word indicator.
    ("NGC 3377", [CAP, CAP, L["n"], L["g"], L["c"], SP, NUM, L["c"], L["c"], L["g"], L["g"]]),
    # Digits end the capitals word; the number sign starts numeric mode.
    ("GRG1", [CAP, CAP, L["g"], L["r"], L["g"], NUM, L["a"]]),
    ("M89", [CAP, L["m"], NUM, L["h"], L["i"]]),
    ("0", [NUM, L["j"]]),
    # A letter a-j straight after a number needs the grade 1 indicator...
    ("1a", [NUM, L["a"], G1, L["a"]]),
    # ...but k-z do not, since they cannot be digits.
    ("1k", [NUM, L["a"], L["k"]]),
    # Decimal point stays in numeric mode.
    ("3.5", [NUM, L["c"], [0, 1, 0, 0, 1, 1], L["e"]]),
    ("Hoag's", [CAP, L["h"], L["o"], L["a"], L["g"], [0, 0, 1, 0, 0, 0], L["s"]]),
    ("Wolf-Rayet", [CAP, L["w"], L["o"], L["l"], L["f"], [0, 0, 1, 0, 0, 1],
                    CAP, L["r"], L["a"], L["y"], L["e"], L["t"]]),
])
def test_ueb_cells(text, expected):
    assert ueb_cells(text) == expected


@pytest.mark.parametrize("text", ["Sgr A*", "Café", "50%", "x²"])
def test_unsupported_characters_raise(text):
    with pytest.raises(ValueError, match="braille"):
        ueb_cells(text)


def test_braille_mesh_dot_count():
    # "SBb" = cap(1 dot) s(3) cap(1) b(2) b(2) = 9 dots, cells 6 mm apart.
    m = core.generate_braille_mesh("SBb")
    assert m.split_bodies().n_blocks == 9
    # First cell (capital sign) is dot 6 only (right column, x = 2.5); last
    # cell "b" is dots 1-2 only (left column, x = 4 * 6); dots are 1.5 mm across.
    assert m.bounds[0] == pytest.approx(2.5 - 0.75, abs=0.05)
    assert m.bounds[1] == pytest.approx(4 * 6.0 + 0.75, abs=0.05)
    assert m.bounds[5] == pytest.approx(core.DEFAULT_BRAILLE_DOT_HEIGHT_MM, abs=0.01)


# --- Labels on the model --------------------------------------------------------

FLAT = np.zeros((50, 50), dtype=np.float32)


def test_braille_too_wide_raises_instead_of_shrinking():
    with pytest.raises(ValueError, match="cannot be scaled down"):
        core.generate_mesh(FLAT, longest_side_mm=40, label_text="Giant Radio Galaxy",
                           label_mode="braille", label_area_mm=15)


def test_text_too_wide_is_scaled_to_fit():
    m = core.generate_mesh(FLAT, longest_side_mm=40, label_text="A rather long label",
                           label_mode="text", label_area_mm=15)
    xs = m.vectors[..., 0]
    assert xs.max() - xs.min() == pytest.approx(40 * 49 / 50, abs=0.01)


def test_braille_dot_height_default_and_override():
    def top(**kw):
        m = core.generate_mesh(FLAT, longest_side_mm=100, base_thickness_mm=2.0,
                               label_text="ab", label_mode="braille", label_area_mm=15, **kw)
        return m.vectors[..., 2].max()
    assert top() == pytest.approx(2.0 + core.DEFAULT_BRAILLE_DOT_HEIGHT_MM, abs=0.01)
    assert top(label_depth_mm=0.6) == pytest.approx(2.6, abs=0.01)


def test_label_errors_propagate():
    with pytest.raises(ValueError):
        core.generate_mesh(FLAT, longest_side_mm=100, label_text="Sgr A*", label_mode="braille")
    with pytest.raises(ValueError, match="label_mode"):
        core.generate_mesh(FLAT, longest_side_mm=100, label_text="x", label_mode="Braille")


# --- Preprocessing ------------------------------------------------------------

def ramp(n=40):
    y, x = np.mgrid[0:n, 0:n]
    return (1.0 + (x + y) / (2.0 * n)).astype(np.float32)


def test_preprocess_does_not_modify_input():
    data = ramp()
    data[5, 5] = np.nan
    before = data.copy()
    core.preprocess_image(data, clip_percentile=1.0, log_scale=True)
    core.preprocess_image(data, clip_percentile=0, downsample_factor=2)
    np.testing.assert_array_equal(data, before)


def test_single_nan_does_not_spread_when_downsampling():
    data = ramp()
    data[20, 20] = np.nan
    out = core.preprocess_image(data, clip_percentile=0, downsample_factor=2)
    # Only the ramp's own minimum corner should sit at 0, not a hole around the NaN.
    assert np.count_nonzero(out == 0) == 1
    assert out[10, 10] > 0.3


def test_blank_region_survives_downsampling():
    data = ramp()
    data[:, :20] = np.nan
    out = core.preprocess_image(data, clip_percentile=0, downsample_factor=2)
    assert out.shape == (20, 20)
    assert np.all(out[:, :9] == 0)
    assert np.all(out[:, 11:] > 0)


def test_all_nan_raises():
    with pytest.raises(ValueError, match="no finite"):
        core.preprocess_image(np.full((10, 10), np.nan, dtype=np.float32))


# --- Mesh -------------------------------------------------------------------------

def test_plaque_is_watertight_and_sized():
    norm = core.preprocess_image(ramp(30), clip_percentile=0)
    m = core.generate_mesh(norm, longest_side_mm=100, max_height_mm=5, base_thickness_mm=2,
                           border_width_mm=5, border_height_mm=1, label_area_mm=10)
    faces = np.hstack([np.full((len(m.vectors), 1), 3), np.arange(len(m.vectors) * 3).reshape(-1, 3)])
    poly = pv.PolyData(m.vectors.reshape(-1, 3).astype(float), faces).clean()
    edges = poly.extract_feature_edges(boundary_edges=True, non_manifold_edges=True,
                                       feature_edges=False, manifold_edges=False)
    assert edges.n_cells == 0
    assert m.vectors[..., 2].min() == 0
    assert m.vectors[..., 2].max() == pytest.approx(7.0, abs=1e-4)


def test_importing_core_does_not_configure_root_logging():
    import importlib
    import logging
    root = logging.getLogger()
    before = list(root.handlers)
    importlib.reload(core)
    assert root.handlers == before


# --- Multi-line labels ------------------------------------------------------------

def _label_dots_z(m, z0):
    """Vertices of the label part of a plaque (everything above the flat strip)."""
    v = m.vectors.reshape(-1, 3)
    return v[v[:, 2] > z0 + 1e-3]


def test_two_braille_lines_are_stacked_in_the_label_area():
    m = core.generate_mesh(FLAT, longest_side_mm=110, base_thickness_mm=2, label_area_mm=24,
                           label_text=["Centaurus A", "Radio"], label_mode="braille")
    dots = _label_dots_z(m, 2.0)
    bottom = m.vectors[..., 1].min()
    # Two lines 10 mm apart, the block centred in the 24 mm strip: line cell
    # centres at 12 +/- 5 mm above the bottom edge, rows 2.5 mm either side.
    assert dots[:, 1].min() - bottom == pytest.approx(7 - 2.5 - 0.75, abs=0.05)
    assert dots[:, 1].max() - bottom == pytest.approx(17 + 2.5 + 0.75, abs=0.05)


def test_label_string_with_newline_equals_list():
    a = core.generate_mesh(FLAT, longest_side_mm=110, label_area_mm=24, label_text="Helix Nebula\nNGC 7293",
                           label_mode="braille")
    b = core.generate_mesh(FLAT, longest_side_mm=110, label_area_mm=24, label_text=["Helix Nebula", "NGC 7293"],
                           label_mode="braille")
    np.testing.assert_allclose(a.vectors, b.vectors)


def test_text_label_stands_full_depth_proud():
    m = core.generate_mesh(FLAT, longest_side_mm=110, base_thickness_mm=2, border_width_mm=3,
                           border_height_mm=1.5, label_area_mm=24, label_text=["Centaurus A", "Visible light"],
                           label_mode="text", label_depth_mm=1.0)
    assert m.vectors[..., 2].max() == pytest.approx(2 + 1.5 + 1.0, abs=1e-4)
    label, _ = core.generate_label_mesh(["Centaurus A", "Visible light"], mode="text", depth_mm=1.0)
    assert label.bounds[4] == pytest.approx(0.0, abs=1e-6)     # stands on the surface...
    assert label.bounds[5] == pytest.approx(1.0, abs=1e-6)     # ...and is fully 1 mm proud


def test_text_lines_share_cap_height_and_baseline():
    # Capitals span exactly the cap height, centred on y=0, whatever else is in the line.
    for text in ["H", "Hg", "HELL"]:     # flat-topped capitals (round ones overshoot slightly)
        m = core.generate_text_mesh(text, cap_height_mm=5.0, depth_mm=1.0)
        assert m.bounds[3] == pytest.approx(2.5, abs=0.01), text
    g = core.generate_text_mesh("Hg", cap_height_mm=5.0)
    assert g.bounds[2] < -2.5 - 1.0       # descender hangs below the cap band
    assert core.generate_text_mesh("ace", 5.0).bounds[3] < 2.5   # x-height letters stay low


def test_text_mesh_is_closed():
    m = core.generate_text_mesh("Radio 3C", cap_height_mm=5.0).clean()
    edges = m.extract_feature_edges(boundary_edges=True, non_manifold_edges=True,
                                    feature_edges=False, manifold_edges=False)
    assert edges.n_cells == 0


def test_wide_text_lines_shrink_together():
    m, _ = core.generate_label_mesh(["A very long first line indeed", "Short"], mode="text",
                                    text_height_mm=5, max_width_mm=60)
    assert m.bounds[1] - m.bounds[0] == pytest.approx(60, abs=0.05)


def test_too_wide_braille_line_is_named():
    with pytest.raises(ValueError, match="'Giant Radio Galaxy 1'"):
        core.generate_label_mesh(["GRG1", "Giant Radio Galaxy 1"], mode="braille", max_width_mm=80)


def test_blank_label_lines_are_ignored():
    m, h = core.generate_label_mesh(["", "  "], mode="braille")
    assert m is None and h == 0


# --- Noise cleaning ------------------------------------------------------------------

from astrotouch import clean  # noqa: E402


def synthetic_radio(seed=1, n=300, noise=1.0):
    rng = np.random.default_rng(seed)
    y, x = np.mgrid[0:n, 0:n]
    img = rng.normal(0, noise, (n, n))
    from scipy.ndimage import gaussian_filter
    img = gaussian_filter(img, 2) / gaussian_filter(rng.normal(0, 1, (n, n)), 2).std() * noise
    img += 200 * np.exp(-((x - 80) ** 2 + (y - 80) ** 2) / (2 * 4 ** 2))      # bright compact source
    img += 2.5 * np.exp(-((x - 200) ** 2 + (y - 190) ** 2) / (2 * 25 ** 2))  # faint 2.5-sigma lobe
    return img.astype(np.float32)


def test_clean_flattens_sky_and_keeps_sources():
    img = synthetic_radio()
    from scipy import ndimage

    def big_islands(im):
        lab, n = ndimage.label(im > 0)
        return int((ndimage.sum(np.ones_like(lab), lab, np.arange(1, n + 1)) >= 20).sum())

    y, x = np.mgrid[0:300, 0:300]
    lobe = (x - 200) ** 2 + (y - 190) ** 2 < 40 ** 2

    # Default: one beam-sized detection scale keeps the compact source, no
    # noise islands, and only a fragment of the 2.5-sigma broad lobe.
    r = clean.clean_image(img, beam_px=(4.7, 4.7))
    out = r.image
    assert out[80, 80] == pytest.approx(200, rel=0.15)   # units of the noise rms
    assert big_islands(out) == 2                     # the source and the lobe fragment
    assert np.mean(out == 0) > 0.9                   # empty sky is exactly flat
    assert np.mean(out[lobe] > 0) < 0.3
    assert r.rms == pytest.approx(1.0, rel=0.15)

    # Adding a broad scale lifts the whole faint lobe above its own noise.
    broad = clean.clean_image(img, beam_px=(4.7, 4.7), smooth_scales_px=(2, 12)).image
    assert big_islands(broad) == 2
    assert np.mean(broad[lobe] > 0) > 0.9
    assert broad[190, 200] > 1.5


def test_clean_handles_blank_regions_and_does_not_modify_input():
    img = synthetic_radio()
    img[:, :40] = np.nan
    before = img.copy()
    r = clean.clean_image(img, beam_px=4.7)
    np.testing.assert_array_equal(img, before)
    assert np.all(r.image[:, :40] == 0)
    assert np.isfinite(r.image).all()


def test_local_noise_suppresses_noisy_half():
    img = synthetic_radio(seed=3)
    img[:, 150:] *= 3                                # noisier right half (and brighter lobe edge)
    img[80, 80] = 0
    glob = clean.clean_image(img, beam_px=4.7)
    loc = clean.clean_image(img, beam_px=4.7, noise_box_px=50)
    from scipy import ndimage
    def islands(im):   # in the noisy half, away from the sharp step at x=150
        lab, n = ndimage.label(im[:, 175:] > 0)
        return int((ndimage.sum(np.ones_like(lab), lab, np.arange(1, n + 1)) >= 20).sum())
    assert islands(glob.image) > islands(loc.image)  # global rms lets noise islands through
    assert islands(loc.image) == 1                   # only the lobe


def test_beam_from_header():
    from astropy.io import fits
    h = fits.Header({"BMAJ": 15.5 / 3600, "BMIN": 11.9 / 3600, "CDELT2": 1.7 / 3600})
    assert clean.beam_from_header(h) == pytest.approx((15.5 / 1.7, 11.9 / 1.7))
    assert clean.beam_from_header(fits.Header({"CDELT2": 1.0})) is None


def test_clean_rejects_noiseless_image():
    with pytest.raises(ValueError, match="noise"):
        clean.clean_image(np.zeros((50, 50)))


# --- Soften (unit-independent stretch) --------------------------------------------

from astrotouch import shaping, sources  # noqa: E402


def test_asinh_without_soften_is_linear_on_small_units_and_soften_fixes_it():
    rng = np.random.default_rng(0)
    y, x = np.mgrid[0:100, 0:100]
    jy = (1e-3 * np.exp(-((x - 50) ** 2 + (y - 50) ** 2) / 200) + rng.normal(0, 1e-5, (100, 100))).astype(np.float32)
    lin = core.preprocess_image(jy, clip_percentile=0)
    plain = core.preprocess_image(jy, clip_percentile=0, asinh_scale=True)
    soft = core.preprocess_image(jy, clip_percentile=0, asinh_scale=True, soften=3)
    assert np.abs(plain - lin).max() < 1e-3           # Jy/beam: asinh does nothing
    assert soft[50, 60] > lin[50, 60] + 0.2           # softened: faint emission lifted
    # Same result whatever the brightness units.
    np.testing.assert_allclose(core.preprocess_image(jy * 1e6, clip_percentile=0, asinh_scale=True, soften=3),
                               soft, atol=1e-5)


# --- Block downsampling and the resolution cap ----------------------------------------

def test_block_reduce_mean_ridge_and_blanks():
    a = np.arange(16, dtype=np.float32).reshape(4, 4)
    np.testing.assert_allclose(shaping.block_reduce(a, 2), [[2.5, 4.5], [10.5, 12.5]])
    np.testing.assert_allclose(shaping.block_reduce(a, 2, ridge=1.0), [[5, 7], [13, 15]])
    a[:2, :2] = np.nan
    a[2, 2] = np.nan
    out = shaping.block_reduce(a, 2)
    assert np.isnan(out[0, 0]) and out[1, 1] == pytest.approx((11 + 14 + 15) / 3)
    assert shaping.block_reduce(np.ones((5, 7)), 2).shape == (2, 3)   # partial blocks cropped


def test_block_mean_averages_noise_away_unlike_linear():
    rng = np.random.default_rng(1)
    noise = rng.normal(0, 1, (400, 400)).astype(np.float32)
    mean = shaping.block_reduce(noise, 8)
    from scipy.ndimage import zoom
    linear = zoom(noise, 1 / 8, order=1)
    assert mean.std() < 0.2 < linear.std()


def test_thin_ridge_survives_with_ridge_mixing():
    img = np.zeros((120, 120), np.float32)
    img[:, 61] = 1.0                                     # a one-pixel filament
    assert shaping.block_reduce(img, 12).max() == pytest.approx(1 / 12, abs=1e-6)
    assert shaping.block_reduce(img, 12, ridge=0.4).max() > 0.4


def test_auto_downsample_factor():
    assert shaping.auto_downsample_factor((5000, 4000), 200, 0.3) == 7     # 0.04 -> 0.28..0.32 mm
    assert shaping.auto_downsample_factor((300, 300), 150, 0.3) == 1
    assert shaping.auto_downsample_factor((5000, 4000), 200, 0) == 1
    assert shaping.auto_downsample_factor((5000, 4000), None, 0.3) == 1


# --- Tactile stretch --------------------------------------------------------------------

def test_tactile_stretch_edge_step_flat_base_and_capped_peak():
    y, x = np.mgrid[0:200, 0:200]
    img = 1000 * np.exp(-np.hypot(x - 100, y - 100) / 15)
    h = shaping.tactile_stretch(img, edge=5.0, contrast=0)
    r = np.hypot(x - 100, y - 100)
    assert np.all(h[img <= 5] == 0)                    # outside the edge: flat base
    assert h[img > 5].min() == pytest.approx(0.2, abs=0.02)   # the outline is a 20% step
    assert np.mean(h == 1) > 0                          # capped peak: a small flat-topped summit
    assert np.mean(h[r < 3] > 0.98) > 0.5


def test_tactile_contrast_boost_lifts_arms():
    y, x = np.mgrid[0:300, 0:300]
    r, t = np.hypot(x - 150, y - 150), np.arctan2(y - 150, x - 150)
    disc = 500 * np.exp(-r / 40) * (1 + 0.5 * np.cos(2 * t - r / 10))
    flat = shaping.tactile_stretch(disc, 2.0, contrast=0)
    boosted = shaping.tactile_stretch(disc, 2.0, contrast=3)
    ring = (r > 60) & (r < 80)
    assert boosted[ring].std() > 1.3 * flat[ring].std()   # arms ~45% more pronounced


def test_tactile_cli_edge_defaults_and_errors():
    img = np.zeros((50, 50), np.float32)
    with pytest.raises(ValueError, match="edge"):
        shaping.tactile_stretch(img, 1.0)
    rng = np.random.default_rng(2)
    y, x = np.mgrid[0:200, 0:200]
    sky = 100 + rng.normal(0, 1, (200, 200)) + 300 * np.exp(-np.hypot(x - 100, y - 100) / 20)
    out = core.preprocess_image(sky.astype(np.float32), clip_percentile=0, tactile=True)
    assert np.mean(out == 0) > 0.5 and out.max() == 1


# --- Compact islands and markers ----------------------------------------------------------

def test_suppress_compact_removes_or_lowers_small_islands():
    h = np.zeros((100, 100), np.float32)
    h[10:12, 10:12] = 0.8          # 2x2 px = 1 x 1 mm at 0.5 mm/px: compact
    h[40:70, 40:70] = 0.6          # 15 mm: real
    out, n = shaping.suppress_compact(h, 0.5, 1.5)
    assert n == 1 and out[10, 10] == 0 and out[50, 50] == pytest.approx(0.6)
    low, _ = shaping.suppress_compact(h, 0.5, 1.5, fraction=0.35, cap=0.15)
    assert low[10, 10] == pytest.approx(0.15)
    same, n = shaping.suppress_compact(np.full((20, 20), 0.5, np.float32), 0.5, 1.5)
    assert n == 0                  # no zero background: refuse rather than guess


def test_marker_dome_sits_on_local_relief():
    h = np.full((100, 100), 0.3, np.float32)
    out = shaping.add_markers(h, [(50, 50)], pixel_mm=0.5, diameter_mm=3, height=0.3)
    assert out[50, 50] == pytest.approx(0.6, abs=1e-6)
    assert out[50, 54] == pytest.approx(0.3)           # 2 mm away: outside the 1.5 mm radius
    with pytest.raises(ValueError, match="outside"):
        shaping.add_markers(h, [(500, 5)], 0.5)


# --- Source removal ------------------------------------------------------------------------

def test_inpaint_fills_without_bumps():
    y, x = np.mgrid[0:100, 0:100]
    smooth = (x + 2 * y).astype(float)
    star = smooth.copy()
    star[45:55, 45:55] += 500
    mask = np.zeros(star.shape, bool)
    mask[43:57, 43:57] = True
    filled = sources.inpaint(star, mask)
    np.testing.assert_allclose(filled, smooth, atol=0.5)    # a plane is harmonic: recovered exactly


def test_remove_points_keeps_protected_nucleus():
    rng = np.random.default_rng(4)
    y, x = np.mgrid[0:200, 0:200]
    img = 50 * np.exp(-np.hypot(x - 100, y - 100) / 30) + rng.normal(0, 0.5, (200, 200))
    for px, py in [(30, 40), (160, 50), (100, 100)]:    # two stars and a sharp nucleus
        img += 60 * np.exp(-((x - px) ** 2 + (y - py) ** 2) / 2)
    out, mask = sources.remove_sources(img, points=True, protect=[(100, 100, 6)])
    assert mask[40, 30] and mask[50, 160] and not mask[100, 100]
    assert out[40, 30] < 10 and out[100, 100] > 100


# --- Picture input -------------------------------------------------------------------------

def test_load_png_flips_to_bottom_up(tmp_path):
    from PIL import Image
    rgb = np.zeros((20, 30, 3), np.uint8)
    rgb[0, 0] = (90, 90, 90)          # top-left in the picture
    Image.fromarray(rgb).save(tmp_path / "pic.png")
    data, header = core.load_image(tmp_path / "pic.png", return_header=True)
    assert header is None and data.shape == (20, 30)
    assert data[19, 0] == 90 and data[0, 0] == 0


# --- Pipeline, painted layers, textures ---------------------------------------------------

from astrotouch import pipeline  # noqa: E402


def galaxy_with_star(n=200):
    rng = np.random.default_rng(5)
    y, x = np.mgrid[0:n, 0:n]
    img = 100 + 400 * np.exp(-np.hypot(x - 100, y - 100) / 20) + rng.normal(0, 1, (n, n))
    img += 3000 * np.exp(-((x - 40) ** 2 + (y - 150) ** 2) / 4)       # a bright star at (40, 150)
    return img.astype(np.float32)


def test_regrid_and_mask_scale():
    assert pipeline.mask_scale((1000, 800)) == 1
    assert pipeline.mask_scale((10527, 10121)) == 6
    assert pipeline.mask_grid_shape((10527, 10121)) == (1755, 1687)
    m = np.zeros((4, 4), bool)
    m[1, 2] = True
    up = pipeline.regrid(m, 3, 1, (12, 12))            # grid cells of 3 raw px -> raw
    assert up[3:6, 6:9].all() and up.sum() == 9
    down = pipeline.regrid(up, 1, 2, (6, 6))           # raw -> 2 px model pixels
    assert down[2, 3] and down.sum() == 2


def test_paint_out_fill_removes_star_and_flatten_zeroes():
    b = pipeline.ModelBuilder(galaxy_with_star())
    s = pipeline.Settings(longest_side_mm=100, min_pixel_mm=0, smooth=0, stretch="tactile", contrast=0)
    paint = np.zeros(b.mask_shape, bool)
    yy, xx = np.indices(paint.shape)
    paint[np.hypot(xx - 40, yy - 150) < 8] = True
    before = b.height_map(s).heights
    filled = b.height_map(s, pipeline.Layers(paint=paint, versions={"paint": 1})).heights
    assert before[150, 40] == pytest.approx(1.0) and filled[150, 40] < 0.1
    flat = b.height_map(pipeline.Settings(**{**s.to_dict(), "paint_mode": "flatten"}),
                        pipeline.Layers(paint=np.ones(b.mask_shape, bool), versions={"paint": 2})).heights
    assert flat.max() == 0


def test_protect_layer_keeps_point_source():
    img = galaxy_with_star()
    b = pipeline.ModelBuilder(img)
    s = pipeline.Settings(longest_side_mm=100, min_pixel_mm=0, smooth=0, remove_points=True, max_point_area=200)
    removed = b.height_map(s).heights
    protect = np.zeros(b.mask_shape, bool)
    protect[140:160, 30:50] = True
    kept = b.height_map(s, pipeline.Layers(protect=protect, versions={"protect": 1})).heights
    assert removed[150, 40] < 0.2 < 0.9 < kept[150, 40]


def test_builder_caches_slow_stages(monkeypatch):
    calls = []
    real = pipeline.sources.remove_sources
    monkeypatch.setattr(pipeline.sources, "remove_sources", lambda *a, **k: calls.append(1) or real(*a, **k))
    b = pipeline.ModelBuilder(galaxy_with_star())
    s = pipeline.Settings(longest_side_mm=100, remove_points=True)
    b.height_map(s)
    b.height_map(pipeline.Settings(**{**s.to_dict(), "max_height_mm": 5, "stretch": "asinh"}))
    assert len(calls) == 1                              # stretch/size changes reuse the removal
    b.height_map(pipeline.Settings(**{**s.to_dict(), "point_sigma": 4}))
    assert len(calls) == 2


def test_texture_patterns():
    for kind in shaping.TEXTURES:
        t = shaping.texture_pattern((200, 200), 0.1, kind, 4.0)
        assert t.min() == 0 and t.max() == pytest.approx(1, abs=1e-3)
    dots = shaping.texture_pattern((400, 400), 0.1, "dots", 4.0)
    from scipy import ndimage
    n = ndimage.label(dots > 0)[1]
    assert 70 < n < 140                                 # ~ (40 mm / 4 mm)^2 x 2/sqrt(3) dots
    h = np.full((100, 100), 0.5, np.float32)
    tm = np.zeros((100, 100), np.uint8)
    tm[:, 50:] = 2
    out = shaping.apply_texture(h, tm, 0.2, 4.0, 0.06)
    assert np.all(out[:, :50] == 0.5) and out[:, 50:].max() == pytest.approx(0.56, abs=1e-3)


def test_texture_in_pipeline_raises_relief_by_texture_height():
    b = pipeline.ModelBuilder(np.zeros((300, 300), np.float32) + np.linspace(0, 1, 300)[None, :])
    s = pipeline.Settings(longest_side_mm=60, min_pixel_mm=0, max_height_mm=10, texture_height_mm=0.6,
                          texture_spacing_mm=4)
    tex = np.zeros(b.mask_shape, np.uint8)
    tex[:, :100] = 3
    base = b.height_map(s).heights
    out = b.height_map(s, pipeline.Layers(texture=tex, versions={"texture": 1})).heights
    diff = (out - base) * 10
    assert diff[:, :100].max() == pytest.approx(0.6, abs=0.01) and np.abs(diff[:, 110:]).max() == 0


def test_layers_round_trip_through_files(tmp_path):
    from astropy.io import fits
    raw_shape = (5000, 3000)
    grid = pipeline.mask_grid_shape(raw_shape)
    tex = np.zeros(grid, np.uint8)
    tex[100:200, 50:80] = 2
    fits.writeto(tmp_path / "t.fits", tex, overwrite=True)
    back = pipeline.load_layer(tmp_path / "t.fits", raw_shape, integer=True)
    np.testing.assert_array_equal(back, tex)
    # A full-resolution mask of a different size is stretched over the image too.
    big = np.zeros(raw_shape, np.uint8)
    big[:2500] = 1
    fits.writeto(tmp_path / "m.fits", big, overwrite=True)
    m = pipeline.load_layer(tmp_path / "m.fits", raw_shape)
    assert m.shape == grid and m[: grid[0] // 2 - 1].all() and not m[grid[0] // 2 + 1:].any()


def test_hillshade_and_plaque_size():
    y, x = np.mgrid[0:100, 0:100]
    dome = np.exp(-((x - 50) ** 2 + (y - 50) ** 2) / 300)
    hs = shaping.hillshade(dome, 0.5, 10)
    assert 0 <= hs.min() and hs.max() <= 1 and hs[50, 50] == pytest.approx(np.sin(np.radians(45)), abs=0.01)
    b = pipeline.ModelBuilder(dome.astype(np.float32))
    s = pipeline.Settings(longest_side_mm=100, min_pixel_mm=0, border_width_mm=5, border_height_mm=1,
                          label_area_mm=20, label=("Dome",), label_mode="braille")
    r = b.height_map(s)
    m = b.mesh(s, r)
    pts = m.vectors.reshape(-1, 3)
    est = pipeline.plaque_size(s, r)
    assert est[0] == pytest.approx(np.ptp(pts[:, 0]), abs=0.6)
    assert est[1] == pytest.approx(np.ptp(pts[:, 1]), abs=0.6)
    assert est[2] == pytest.approx(pts[:, 2].max(), abs=0.05)
    assert pipeline.fits_bed((210, 160), (220, 220)) and pipeline.fits_bed((230, 150), (160, 240))
    assert not pipeline.fits_bed((230, 230), (220, 220))


# --- Gaia service handling (offline: the network is mocked) ---------------------------------

GAIA_CSV = (b"ra,dec,phot_g_mean_mag,parallax_over_error,pmra,pmdec,pmra_error,pmdec_error,ruwe\n"
            b"10.0,20.0,12.0,10.0,5.0,5.0,0.1,0.1,1.0\n"       # a nearby star
            b"10.01,20.01,19.0,0.1,0.01,0.01,0.5,0.5,1.0\n")   # a distant non-star source


class _Resp:
    def __init__(self, data):
        self.data = data

    def read(self):
        return self.data

    def __enter__(self):
        return self

    def __exit__(self, *a):
        return False


def _mock_services(monkeypatch, responses):
    """responses: one per call, either bytes (a reply) or an exception to raise."""
    import urllib.request
    calls = []

    def fake(req, timeout=None):
        calls.append(req.full_url)
        r = responses[len(calls) - 1]
        if isinstance(r, Exception):
            raise r
        return _Resp(r)
    monkeypatch.setattr(urllib.request, "urlopen", fake)
    monkeypatch.setattr(sources, "_FAILED", {})
    return calls


def _http_error(code, msg="Server overloaded"):
    import io
    import urllib.error
    return urllib.error.HTTPError("u", code, "x", {}, io.BytesIO(f"<p>Message: {msg}</p>".encode()))


def test_gaia_falls_back_to_mirror_and_caches(monkeypatch, tmp_path):
    calls = _mock_services(monkeypatch, [_http_error(503), _http_error(503), GAIA_CSV])
    monkeypatch.setattr("time.sleep", lambda s: None)
    stars = sources.gaia_stars(10, 20, 0.1, cache=tmp_path / "g.csv")
    assert len(calls) == 3 and "esa" in calls[0] and "uni-heidelberg" in calls[2]
    assert len(stars) == 1 and stars["phot_g_mean_mag"][0] == 12.0
    assert (tmp_path / "g.csv").read_bytes() == GAIA_CSV


def test_gaia_bad_response_is_not_cached_and_error_is_clear(monkeypatch, tmp_path):
    _mock_services(monkeypatch, [b"<html>maintenance</html>", _http_error(503), _http_error(503)])
    monkeypatch.setattr("time.sleep", lambda s: None)
    with pytest.raises(RuntimeError) as e:
        sources.gaia_stars(10, 20, 0.1, cache=tmp_path / "g.csv")
    msg = str(e.value)
    assert "not a Gaia table" in msg and "HTTP 503 (Server overloaded" in msg and "--gaia-cache" in msg
    assert not (tmp_path / "g.csv").exists()


def test_gaia_failure_is_remembered_briefly(monkeypatch):
    calls = _mock_services(monkeypatch, [_http_error(404), _http_error(404)])
    with pytest.raises(RuntimeError):
        sources.gaia_stars(10, 20, 0.1)
    n = len(calls)
    assert n == 2                                     # 404 is not retried
    with pytest.raises(RuntimeError):
        sources.gaia_stars(10, 20, 0.1)               # within a minute: no new network calls
    assert len(calls) == n


def test_remove_sources_reports_what_was_removed_and_why():
    rng = np.random.default_rng(6)
    y, x = np.mgrid[0:200, 0:200]
    img = 50 * np.exp(-np.hypot(x - 100, y - 100) / 30) + rng.normal(0, 0.5, (200, 200))
    for px, py in [(30, 40), (160, 50)]:
        img += 60 * np.exp(-((x - px) ** 2 + (y - py) ** 2) / 2)
    paint = np.zeros(img.shape, bool)
    paint[150:160, 150:160] = True
    out, mask, parts = sources.remove_sources(img, points=True, extra_masks=[(20, 180, 5)], extra_mask=paint,
                                              protect=[(100, 100, 6)], parts=True)   # keep the nucleus
    assert parts["n_points"] == 2 and parts["points"][40, 30] and parts["points"][50, 160]
    assert parts["discs"][180, 20] and parts["painted"][155, 155] and parts["n_stars"] == 0
    assert not (parts["points"] & (parts["discs"] | parts["painted"])).any()   # each pixel has one reason
    np.testing.assert_array_equal(mask, parts["points"] | parts["discs"] | parts["painted"])


def test_pipeline_preview_carries_removed_outlines():
    b = pipeline.ModelBuilder(galaxy_with_star(400))
    s = pipeline.Settings(longest_side_mm=100, min_pixel_mm=0.5, remove_points=True, max_point_area=200)
    r = b.height_map(s, with_display=True)
    assert r.factor == 2 and r.removed["points"].shape == r.heights.shape
    assert r.removed["n_points"] >= 1 and r.removed["points"][75, 20]          # the star at raw (40, 150)
    assert b.height_map(s).removed is None                                      # only for previews


# --- Shapes and labels on the image ----------------------------------------------------------

def test_raised_line_has_requested_width_and_height():
    h = np.zeros((100, 200), np.float32)
    out = shaping.add_shapes(h, [dict(type="line", points=[(20, 50), (180, 50)], width_mm=1.2, height=0.1)], 0.2)
    assert out.max() == pytest.approx(0.1, abs=1e-6)
    col = out[:, 100]
    assert np.count_nonzero(col) in (5, 6, 7)          # 1.2 mm at 0.2 mm/px ~ 6 px
    assert out[50, 10] == 0 and out[50, 190] == 0      # stops at the end points


def test_arrow_head_and_ring():
    h = np.zeros((200, 200), np.float32)
    out = shaping.add_shapes(h, [dict(type="line", points=[(20, 100), (180, 100)], arrow=True,
                                      width_mm=1.0, height=0.1)], 0.2)
    head = max(np.count_nonzero(out[:, c]) for c in range(155, 181))
    assert head >= 2.5 * np.count_nonzero(out[:, 60])  # arrowhead clearly wider than the shaft
    assert out[100, 179] > 0 and out[100, 182] == 0     # its tip ends at the end point
    ring = shaping.add_shapes(h, [dict(type="circle", x=100, y=100, radius_mm=10, width_mm=1.0, height=0.1)], 0.2)
    assert ring[100, 100] == 0 and ring[100, 150] == pytest.approx(0.1)        # hollow, ring at r = 50 px
    disc = shaping.add_shapes(h, [dict(type="circle", x=100, y=100, radius_mm=10, filled=True, height=0.1)], 0.2)
    assert disc[100, 100] == pytest.approx(0.1)


def test_shapes_follow_relief_and_do_not_stack():
    h = np.full((100, 100), 0.4, np.float32)
    out = shaping.add_shapes(h, [dict(type="line", points=[(10, 50), (90, 50)], height=0.1),
                                 dict(type="line", points=[(50, 10), (50, 90)], height=0.1)], 0.2)
    assert out.max() == pytest.approx(0.5, abs=1e-6)  # crossing lines: 0.4 + 0.1, not + 0.2


def _flat_builder():
    y, x = np.mgrid[0:400, 0:400]
    return pipeline.ModelBuilder((np.exp(-np.hypot(x - 200, y - 200) / 60)).astype(np.float32))


def test_image_label_sits_on_flat_pad_at_the_right_place():
    b = _flat_builder()
    s = pipeline.Settings(longest_side_mm=100, min_pixel_mm=0, smooth=0, base_thickness_mm=2, max_height_mm=10,
                          border_width_mm=4, border_height_mm=1, label_area_mm=15,
                          annotations=({"type": "label", "x": 100, "y": 300, "text": "Core"},))
    r = b.height_map(s)
    (m, gx, gy, level, (c0, r0, c1, r1), text) = r.image_labels[0]
    assert text == "Core" and level == 0
    assert np.all(r.heights[r0:r1, c0:c1] == 0)        # flattened pad
    import dataclasses
    plain = len(b.mesh(s, dataclasses.replace(r, image_labels=[])).vectors)
    dots = b.mesh(s, r).vectors[plain:].reshape(-1, 3)   # the label's triangles come last
    # Dots stand on the base (2 mm) up to 2.8 mm, centred where pixel (100, 300) lands on the model:
    # 0.25 mm pixels, 16 border px each side, a 60 px label strip below.
    assert dots[:, 2].min() == pytest.approx(2.0, abs=1e-4)
    assert dots[:, 2].max() == pytest.approx(2.8, abs=0.01)
    cx = (100 + 16 - (400 + 32 - 1) / 2) * 0.25
    cy = (300 + 16 + 60 - (400 + 32 + 60 - 1) / 2) * 0.25
    assert (dots[:, 0].min() + dots[:, 0].max()) / 2 == pytest.approx(cx, abs=1.5)   # dots, not cells, so ~
    assert (dots[:, 1].min() + dots[:, 1].max()) / 2 == pytest.approx(cy, abs=1.5)


def test_raised_pad_and_label_off_image():
    b = _flat_builder()
    s = pipeline.Settings(longest_side_mm=100, min_pixel_mm=0, smooth=0, image_label_pad="raised",
                          annotations=({"type": "label", "x": 200, "y": 200, "text": "Peak"},))
    r = b.height_map(s)
    _, _, _, level, (c0, r0, c1, r1), _ = r.image_labels[0]
    assert level == pytest.approx(1.0) and np.all(r.heights[r0:r1, c0:c1] == level)
    off = pipeline.Settings(longest_side_mm=100, annotations=({"type": "label", "x": 5, "y": 200, "text": "Edge"},))
    with pytest.raises(ValueError, match="runs off the image"):
        b.height_map(off)


def test_leader_line_stops_at_label_pad():
    b = _flat_builder()
    s = pipeline.Settings(longest_side_mm=100, min_pixel_mm=0, smooth=0, image_label_mode="text",
                          annotations=({"type": "line", "points": [[100, 300], [200, 200]], "arrow": True},
                                       {"type": "label", "x": 100, "y": 300, "text": "Core"}))
    r = b.height_map(s)
    _, _, _, _, (c0, r0, c1, r1), _ = r.image_labels[0]
    assert np.all(r.heights[r0:r1, c0:c1] == 0)
    assert r.heights[r1 + 5, c1 + 5] > 0               # the line continues outside the pad


def test_cli_reads_annotations_json(tmp_path):
    import json
    import subprocess
    from astropy.io import fits
    y, x = np.mgrid[0:200, 0:200]
    fits.writeto(tmp_path / "in.fits", np.exp(-np.hypot(x - 100, y - 100) / 30).astype(np.float32))
    ann = [{"type": "circle", "x": 100, "y": 100, "radius_mm": 15},
           {"type": "label", "x": 100, "y": 30, "text": "Ring"}]
    (tmp_path / "a.json").write_text(json.dumps(ann))
    gui_json = {"settings": {"annotations": ann}}
    (tmp_path / "g.astrotouch.json").write_text(json.dumps(gui_json))
    root = Path(__file__).resolve().parents[1]
    for src in ("a.json", "g.astrotouch.json"):
        r = subprocess.run([sys.executable, str(root / "cli.py"), str(tmp_path / "in.fits"), str(tmp_path / "o.stl"),
                            "--longest-side", "80", "--annotations", str(tmp_path / src)],
                           capture_output=True, text=True, cwd=root)
        assert r.returncode == 0, r.stderr
        assert "2 shape(s) and label(s)" in r.stderr


def test_label_pad_slopes_into_relief_and_warns():
    h = np.zeros((300, 300), np.float32)
    h[100:200, 160:260] = 0.8                          # a flat-topped source right of the label
    b = pipeline.ModelBuilder(h)
    s = pipeline.Settings(longest_side_mm=75, min_pixel_mm=0, smooth=0, stretch="linear", clip=0,
                          annotations=({"type": "label", "x": 150, "y": 150, "text": "Ab"},))
    r = b.height_map(s)
    _, _, _, _, (c0, r0, c1, r1), _ = r.image_labels[0]
    row = r.heights[150]
    assert np.all(row[c0:c1] == 0)                     # flat pad
    steps = np.diff(row[c1 - 1:c1 + 8])
    assert steps.max() < 0.5                           # a slope back up to 0.8, not a 0.8 cliff
    assert any("covers relief" in w for w in r.warnings)
    clear = b.height_map(pipeline.Settings(**{**s.to_dict(), "annotations": (
        {"type": "label", "x": 60, "y": 260, "text": "Ab"},)}))
    assert clear.warnings == []


def test_warns_when_lines_are_too_narrow_for_the_pixels():
    b = pipeline.ModelBuilder(np.zeros((100, 100), np.float32) + np.eye(100, dtype=np.float32))
    ann = ({"type": "line", "points": [[10, 10], [90, 90]]},)
    coarse = b.height_map(pipeline.Settings(longest_side_mm=60, min_pixel_mm=0, annotations=ann))   # 0.6 mm px
    fine = b.height_map(pipeline.Settings(longest_side_mm=20, min_pixel_mm=0, annotations=ann))     # 0.2 mm px
    assert any("jagged" in w for w in coarse.warnings) and not fine.warnings


# --- Templates and the stretch chart ---------------------------------------------------------------

def test_templates_apply_their_own_fields_and_reset_them():
    s = pipeline.Settings(longest_side_mm=77, label=("Keep me",), mask_stars=True, contrast=9)
    radio = pipeline.apply_template(s, recipe="radio")
    assert radio.clean and radio.stretch == "asinh"
    assert not radio.mask_stars and radio.contrast == 1.5          # recipe fields reset to defaults
    assert radio.longest_side_mm == 77 and radio.label == ("Keep me",)   # plaque fields and labels untouched
    both = pipeline.apply_template(radio, plaque="small")
    assert both.longest_side_mm == 90 and both.label_area_mm == 12 and both.clean
    with pytest.raises(ValueError, match="Unknown template"):
        pipeline.apply_template(s, recipe="nope")


def test_every_template_builds_a_model():
    b = pipeline.ModelBuilder(galaxy_with_star(300))
    for plaque in pipeline.PLAQUE_FORMATS:
        for recipe in pipeline.RECIPES:
            s = pipeline.apply_template(pipeline.Settings(), plaque=plaque, recipe=recipe)
            s = pipeline.Settings(**{**s.to_dict(), "mask_stars": False, "edge": 150.0 if s.stretch == "tactile" else None})
            r = b.height_map(s)
            assert np.isfinite(r.heights).all() and r.heights.max() > 0, (plaque, recipe)


def test_saved_template_round_trip_excludes_per_image_content():
    s = pipeline.Settings(stretch="tactile", contrast=4, markers=((1, 2),), label=("M74",),
                          annotations=({"type": "label", "x": 1, "y": 2, "text": "x"},))
    d = pipeline.template_dict(s)
    assert "markers" not in d and "annotations" not in d and "label" not in d
    back = pipeline.apply_template(pipeline.Settings(), custom=d)
    assert back.contrast == 4 and back.stretch == "tactile" and back.markers == ()


def test_stretch_curve_reports_flat_base():
    b = pipeline.ModelBuilder(galaxy_with_star(300))
    r = b.height_map(pipeline.Settings(stretch="tactile", edge=150.0, longest_side_mm=100), with_display=True)
    c = r.curve
    assert 0.5 < c["flat_fraction"] < 1 and len(c["rank"]) == len(c["height"]) <= 3000
    faint = c["height"][c["rank"] < 50]
    assert np.all(faint == 0)                                       # the fainter half is flat sky


def test_cli_templates_and_overrides(tmp_path):
    import subprocess
    from astropy.io import fits
    y, x = np.mgrid[0:200, 0:200]
    fits.writeto(tmp_path / "in.fits", np.exp(-np.hypot(x - 100, y - 100) / 30).astype(np.float32))
    root = Path(__file__).resolve().parents[1]
    r = subprocess.run([sys.executable, str(root / "cli.py"), str(tmp_path / "in.fits"), str(tmp_path / "o.stl"),
                        "--plaque", "small", "--max-height", "7"], capture_output=True, text=True, cwd=root)
    assert r.returncode == 0, r.stderr
    assert "Using plaque 'small'; overridden by max_height_mm" in r.stderr
    assert "Model size: 96 x" in r.stderr                            # 90 mm + 3 mm border each side
    listing = subprocess.run([sys.executable, str(root / "cli.py"), "--list-templates"], capture_output=True,
                             text=True, cwd=root)
    assert "press_picture" in listing.stdout and "standard" in listing.stdout
