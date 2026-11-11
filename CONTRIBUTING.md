# Contributing to AstroTouch

Thank you for helping make astronomy something everyone can feel. You don't need to write code to contribute.

## Ways to help

* **Try the models with people who are blind or have low vision**, and tell us what could and couldn't be felt.
  Use the [feedback form](../../issues/new?template=accessibility_feedback.yml). This is the most useful contribution of all.
* **Check braille.** Labels are machine-encoded Unified English Braille (uncontracted). A braille reader's check of
  the encoder (`core.ueb_cells`) or of a printed label is very welcome.
* **Share models and data.** A model you made (attach the `.astrotouch.json` the GUI saves) or an object and dataset
  that would make a good one: use [this form](../../issues/new?template=model_or_data.yml).
* **Report bugs and suggest features** with the issue forms.
* **Improve the code, the help text or the docs.** See below.

## Development setup

```bash
git clone https://github.com/enzoperesafonso/AstroTouch.git
cd AstroTouch
python -m venv .venv && source .venv/bin/activate    # Windows: .venv\Scripts\activate
pip install -e ".[gui,dev]"
pytest                                               # the tests
ruff check src tests cli.py gui.py                   # the linter
astrotouch-gui                                       # the desktop app (or: python gui.py)
```

Python 3.10 or later. The tests use small synthetic images, so no data downloads are needed; the GUI test is
skipped when the GUI extras aren't installed. Continuous integration runs the tests and the linter on Linux,
macOS and Windows for every pull request.

## How the code fits together

Everything goes through one pipeline (`src/astrotouch/pipeline.py`), shared by the GUI, the command line and the
Python API:

```
load -> remove stars / painted areas -> clean noise -> downsample + stretch
     -> compact islands -> flatten paint-out -> smooth -> markers -> textures
     -> raised shapes -> label pads -> mesh (border, label strip, labels)
```

| Module | What it does |
| :--- | :--- |
| `pipeline.py` | `Settings` (every option, with the shared defaults), `ModelBuilder` (runs the pipeline and caches slow stages), templates, painted layers |
| `core.py` | loading images, preprocessing and stretches, mesh generation, braille and text labels, the command line |
| `clean.py` | noise cleaning for radio maps |
| `sources.py` | Gaia star removal, compact-source removal, harmonic inpainting |
| `shaping.py` | block downsampling, the tactile stretch, compact islands, markers, textures, raised shapes, hillshade |
| `gui.py` | the desktop GUI (PySide6): settings sidebar, 2D paint view, 3D view |
| `help_text.py` | the GUI's help panels, getting-started guide and tool hints |

### Adding a setting

A new option belongs in all three front ends:

1. a field in `pipeline.Settings` (with its default), used in `ModelBuilder`;
2. a command-line flag in `core.main()`, passed into `Settings`, and added to the template override table there;
3. a widget in `gui.py`: create it in `_build_sidebar`, read it in `settings()`, and set it in `apply_settings()`;
4. a line of help in `help_text.py`, a row in the README's reference tables, and a test.

### Adding a template

Templates are plain dictionaries in `pipeline.py`: `PLAQUE_FORMATS` (size, border, label strip) and `RECIPES`
(how the image is prepared). Give each a `short` name, a `title` and an `about` sentence, which appear in the GUI and
`--list-templates`. `test_every_template_builds_a_model` checks that it works.

## Guidelines

* **Design for touch.** Changes to the models should follow tactile-graphics guidance: features about 0.5-1 mm
  high and about 3 mm apart, edges and steps rather than gentle slopes, nothing that could be mistaken for braille,
  and braille never scaled. Say in the PR how the change will feel.
* **Be honest about what isn't data.** Markers, textures, shapes and exaggerations should be documented as such.
* **Data and images.** Credit every image and dataset and follow its licence (for example, MeerKAT data are
  CC BY-NC). Don't commit large data files; tests use synthetic images.
* **Tests.** Add a test for new behaviour, and for bugs, a test that fails without the fix.
* **Keep pull requests focused**, one change per PR, with a short description of what and why.

## Code of conduct

This project follows the [Code of Conduct](CODE_OF_CONDUCT.md). Please be kind; many contributors and users are
new to astronomy, to programming, or to tactile design.
