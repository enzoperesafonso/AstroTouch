# Changelog

## 0.2.0 (unreleased)

### Added
* **Desktop GUI rewrite:** every option in collapsible sections with plain-language help; a live 2D preview
  (shaded relief or input image) rebuilt in the background; plaque size and print-bed check; painting tools
  (paint-out, protect, texture, erase, markers); labels, lines, arrows and rings on the image; outlines of removed
  stars and sources; a chart of what the stretch does; a getting-started guide and tool hints; drag and drop.
  Saving writes the painted layers and all settings beside the STL.
* **Templates:** plaque formats (large, standard, small) and image recipes (radio map, radio mosaic, smooth or
  knotty galaxy, nebula, press picture), in the GUI and as `--plaque`, `--recipe` and `--template`.
* **Noise cleaning** for radio maps (`--clean`), with global or local noise and multi-scale detection.
* **Star and source removal:** Gaia DR3 foreground stars (with a fallback mirror, retries and caching),
  compact sources, hand-placed discs and painted areas, filled by harmonic inpainting.
* **Tactile stretch** (`--tactile`): flat sky, an outline step, a capped peak and a local-contrast boost.
* `--soften` for log and asinh stretches of data with small values (e.g. Jy/beam).
* Block-mean downsampling with ridge preservation, and an automatic resolution cap (`--min-pixel`).
* Compact-island removal, marker domes, tactile textures, raised shapes and labels on the image.
* Multi-line labels; picture input (PNG, JPEG, TIFF); `--bed` size check.
* A shared pipeline (`astrotouch.pipeline`) for the GUI, the command line and Python; packaging with `astrotouch`
  and `astrotouch-gui` commands; tests and continuous integration.

### Changed
* Braille dots default to 0.8 mm high (standard braille), raised text to 1.5 mm.
* `--longest-side` defaults to 150 mm and `--smooth` to 1.0; downsampling uses a block mean by default
  (`--downsample-method linear` for the old behaviour).
* Python 3.10 or later.

### Fixed
* Braille: mixed-case words ("SBb") were encoded as all capitals; numbers followed by letters a-j lacked the grade 1
  indicator; unsupported characters were silently dropped; too-wide braille was shrunk.
* Raised text stood only half its depth above the surface, and lines with descenders were misaligned.
* Label errors were swallowed, saving models without their label.
* NaNs spread when downsampling; the input array was modified in place; the GUI preview was upside down.
* Importing the library reconfigured the caller's logging.
