"""Smoke test of the desktop GUI, offscreen. Skipped where the GUI extras aren't installed."""
import os
import sys
import time

import numpy as np
import pytest

pytest.importorskip("PySide6")
pytest.importorskip("pyvistaqt")
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")


@pytest.fixture(scope="module")
def app(tmp_path_factory):
    from PySide6 import QtWidgets, QtCore
    # Keep the tests out of the user's real preferences, and skip the first-run guide.
    prefs = tmp_path_factory.mktemp("prefs") / "astrotouch.ini"
    os.environ["ASTROTOUCH_SETTINGS_FILE"] = str(prefs)
    QtCore.QSettings(str(prefs), QtCore.QSettings.IniFormat).setValue("seen_intro", True)
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
    yield app
    os.environ.pop("ASTROTOUCH_SETTINGS_FILE", None)


def _wait(app, w, limit=60):
    t = time.time()
    while time.time() - t < limit:
        app.processEvents()
        busy = w.timer.isActive() or (w.preview_worker is not None and w.preview_worker.isRunning())
        if not busy and time.time() - t > 0.5:
            break
        time.sleep(0.05)
    app.processEvents()


def test_open_preview_template_and_annotations(app, tmp_path):
    from astropy.io import fits
    from astrotouch.gui import AstroTouchWindow
    y, x = np.mgrid[0:200, 0:200]
    fits.writeto(tmp_path / "g.fits", (100 + 500 * np.exp(-np.hypot(x - 100, y - 100) / 25)).astype(np.float32))
    w = AstroTouchWindow()
    w.load_file(str(tmp_path / "g.fits"))
    w.use_template(plaque="small", recipe="plain")
    w.canvas.annotations = [{"type": "circle", "x": 100.0, "y": 100.0, "radius_mm": 8}]
    w.schedule_preview()
    _wait(app, w)
    assert w.last_result is not None, w.status.text()
    result, settings = w.last_result
    assert settings.longest_side_mm == 90 and result.curve is not None
    assert "Plaque" in w.lbl_summary.text()
    w.close()
