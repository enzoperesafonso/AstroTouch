#!/usr/bin/env python3
"""
AstroTouch desktop GUI.

Left: every model setting, grouped. Right: a 2D view of the height map (with
brushes to paint out, protect and texture areas and to place markers) and the
3D model. The 2D preview is rebuilt in the background a moment after any
change; "Process 3D model" builds the mesh.
"""
import hashlib
import json
import os
import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np
from astropy.io import fits
from PySide6 import QtWidgets, QtCore, QtGui
import pyvista as pv
from pyvistaqt import QtInteractor

try:
    from . import pipeline, shaping, help_text
    from .core import DEFAULT_TEXT_DEPTH_MM, DEFAULT_BRAILLE_DOT_HEIGHT_MM, IMAGE_SUFFIXES
except ImportError:
    import pipeline
    import shaping
    import help_text
    from core import DEFAULT_TEXT_DEPTH_MM, DEFAULT_BRAILLE_DOT_HEIGHT_MM, IMAGE_SUFFIXES

OPENABLE = tuple(IMAGE_SUFFIXES) + (".fits", ".fit", ".fts", ".fz", ".gz")

DEFAULTS = pipeline.Settings()
PREVIEW_DELAY_MS = 350
PAINT_RGBA = (255, 50, 50, 120)
PROTECT_RGBA = (40, 220, 90, 120)
TEXTURE_RGBA = {1: (255, 170, 0, 130), 2: (0, 190, 255, 130), 3: (190, 90, 255, 130)}
# Outlines of what source removal took out automatically, by reason.
REMOVED_RGBA = {"stars": (0, 230, 255, 255), "points": (255, 70, 220, 255), "discs": (255, 255, 255, 255)}
TOOLS = [("pan", "Pan", "Drag to move; scroll to zoom (any tool)"),
         ("paint", "Paint out", "Paint areas to remove: filled from their surroundings, or flattened"),
         ("protect", "Protect", "Paint areas that star and source removal must never touch"),
         ("texture", "Texture", "Paint a tactile texture onto the model"),
         ("erase", "Erase", "Erase paint-out, protect and texture"),
         ("marker", "Marker", "Click to place a marker dome; right-click one to remove it"),
         ("label", "Label", "Click to place a braille or text label on the image; double-click one to edit it"),
         ("line", "Line", "Click points of a raised line; double-click or Enter to finish, Esc to cancel"),
         ("arrow", "Arrow", "Like Line, with an arrowhead at the last point"),
         ("circle", "Circle", "Drag from the centre to draw a raised ring")]
ANNOTATION_TOOLS = ("label", "line", "arrow", "circle")
SHAPE_COLOUR = QtGui.QColor("#ffd24d")
LABEL_COLOUR = QtGui.QColor("#9be7ff")


# --- Background work ------------------------------------------------------------------

class PreviewWorker(QtCore.QThread):
    done = QtCore.Signal(object, object)       # (HeightResult, Settings)
    failed = QtCore.Signal(str)

    def __init__(self, builder, settings, layers):
        super().__init__()
        self.builder, self.settings, self.layers = builder, settings, layers

    def run(self):
        try:
            self.done.emit(self.builder.height_map(self.settings, self.layers, with_display=True), self.settings)
        except Exception as e:  # shown in the status line
            self.failed.emit(str(e))


class BuildWorker(QtCore.QThread):
    done = QtCore.Signal(object)
    failed = QtCore.Signal(str)
    progress = QtCore.Signal(str)

    def __init__(self, builder, settings, layers):
        super().__init__()
        self.builder, self.settings, self.layers = builder, settings, layers

    def run(self):
        try:
            self.progress.emit("Building height map...")
            result = self.builder.height_map(self.settings, self.layers)
            self.progress.emit("Building mesh...")
            stl_mesh = self.builder.mesh(self.settings, result)
            path = Path(tempfile.gettempdir()) / f"astrotouch_preview_{os.getpid()}.stl"
            stl_mesh.save(str(path))
            pts = stl_mesh.vectors.reshape(-1, 3)
            self.done.emit({'path': path, 'faces': len(stl_mesh.vectors),
                            'size': (float(np.ptp(pts[:, 0])), float(np.ptp(pts[:, 1])), float(pts[:, 2].max()))})
        except Exception as e:
            self.failed.emit(str(e))


# --- 2D view and painting ----------------------------------------------------------------

class PaintCanvas(QtWidgets.QWidget):
    """
    Shows the height map (or the input) and lets the user paint layers on the
    mask grid. Display pixels are height-map pixels, drawn top row first; the
    arrays underneath keep FITS orientation (row 0 at the bottom).
    """
    layersChanged = QtCore.Signal()
    hint = QtCore.Signal(str)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setMouseTracking(True)
        self.setMinimumSize(320, 320)
        self.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.base = None
        self.overlay = None
        self.disp_shape = None
        self.factor = 1
        self.pixel_mm = 1.0
        self.raw_shape = None
        self.ms = 1
        self.paint = self.protect = self.texture = None
        self.markers = []
        self.annotations = []      # dicts in raw pixels (see pipeline.Settings.annotations)
        self.label_boxes = []      # ((c0, r0, c1, r1), text) on the height-map grid, from the last preview
        self.shape_width_mm = DEFAULTS.shape_width_mm
        self.fill_circles = False
        self._pending = []         # points of a line being drawn
        self._circle = None        # centre of a circle being drawn
        self.removed = None
        self.show_removed = True
        self.versions = {"paint": 0, "protect": 0, "texture": 0}
        self.tool = "paint"
        self.brush_mm = 3.0
        self.texture_kind = 1
        self.marker_size_mm = DEFAULTS.marker_size_mm
        self.zoom = 1.0
        self.pan = QtCore.QPointF(0, 0)
        self._undo = []
        self._last = None
        self._drag = None
        self._mouse = None
        self._changed = set()

    # Layers -----------------------------------------------------------------------
    def reset_layers(self, raw_shape):
        self.raw_shape = tuple(raw_shape)
        self.ms = pipeline.mask_scale(raw_shape)
        shape = pipeline.mask_grid_shape(raw_shape)
        self.paint = np.zeros(shape, bool)
        self.protect = np.zeros(shape, bool)
        self.texture = np.zeros(shape, np.uint8)
        self.markers = []
        self.annotations = []
        self.label_boxes = []
        self._pending = []
        self._undo = []
        for k in self.versions:
            self.versions[k] += 1
        self._update_overlay()

    def layers(self):
        """A snapshot of the layers for a background build."""
        if self.paint is None:
            return None
        return pipeline.Layers(paint=self.paint.copy(), protect=self.protect.copy(),
                               texture=self.texture.copy(), versions=dict(self.versions))

    def has_layers(self):
        return self.paint is not None and (self.paint.any() or self.protect.any() or self.texture.any())

    def clear(self, what):
        self._push_undo()
        if what == "markers":
            self.markers = []
        elif what == "annotations":
            self.annotations = []
            self._pending = []
        else:
            arr = getattr(self, what)
            arr[...] = 0
            self.versions[what] += 1
        self._update_overlay()
        self.layersChanged.emit()

    def _push_undo(self):
        if self.paint is None:
            return
        import copy
        self._undo.append((self.paint.copy(), self.protect.copy(), self.texture.copy(), list(self.markers),
                           copy.deepcopy(self.annotations)))
        del self._undo[:-15]

    def undo(self):
        if not self._undo:
            return
        self.paint, self.protect, self.texture, self.markers, self.annotations = self._undo.pop()
        self._pending = []
        for k in self.versions:
            self.versions[k] += 1
        self._update_overlay()
        self.layersChanged.emit()

    # Image ----------------------------------------------------------------------------
    def set_removed(self, removed):
        """Masks (height-map grid, FITS orientation) of what source removal took out."""
        self.removed = removed
        self._update_overlay()

    def set_show_removed(self, on):
        self.show_removed = on
        self._update_overlay()

    def set_image(self, gray, factor, pixel_mm):
        """gray: display-oriented (top row first) array in [0, 1]."""
        g = np.ascontiguousarray((np.clip(gray, 0, 1) * 255).astype(np.uint8))
        h, w = g.shape
        refit = self.disp_shape != (h, w)
        self.base = QtGui.QImage(g.data, w, h, w, QtGui.QImage.Format_Grayscale8).copy()
        self.disp_shape = (h, w)
        self.factor, self.pixel_mm = factor, pixel_mm
        if refit:
            self.zoom, self.pan = 1.0, QtCore.QPointF(0, 0)
        self._update_overlay()

    def _update_overlay(self):
        if self.disp_shape is None or self.paint is None:
            self.overlay = None
            self.update()
            return
        rgba = np.zeros(self.disp_shape + (4,), np.uint8)
        for layer, colour in ((self.protect, PROTECT_RGBA), (self.paint, PAINT_RGBA)):
            if layer.any():
                rgba[np.flipud(pipeline.regrid(layer, self.ms, self.factor, self.disp_shape))] = colour
        if self.texture.any():
            t = np.flipud(pipeline.regrid(self.texture, self.ms, self.factor, self.disp_shape))
            for k, colour in TEXTURE_RGBA.items():
                rgba[t == k] = colour
        if self.show_removed and self.removed is not None:
            from scipy import ndimage
            for key, colour in REMOVED_RGBA.items():
                m = self.removed.get(key)
                if m is not None and m.shape == self.disp_shape and m.any():
                    # A ring just outside each removed area, so in "Input image"
                    # mode the removed star itself stays visible inside it.
                    ring = ndimage.binary_dilation(m) & ~m
                    rgba[np.flipud(ring)] = colour
        h, w = self.disp_shape
        self.overlay = QtGui.QImage(rgba.data, w, h, 4 * w, QtGui.QImage.Format_RGBA8888).copy()
        self.update()

    # Coordinates --------------------------------------------------------------------
    def _scale(self):
        h, w = self.disp_shape
        return min(self.width() / w, self.height() / h) * self.zoom

    def _origin(self):
        h, w = self.disp_shape
        s = self._scale()
        return QtCore.QPointF((self.width() - w * s) / 2, (self.height() - h * s) / 2) + self.pan

    def _to_disp(self, pt):
        o, s = self._origin(), self._scale()
        return (pt.x() - o.x()) / s, (pt.y() - o.y()) / s

    def _disp_to_raw(self, dx, dy):
        """Continuous display coords (y down) to raw pixel-centre coords (FITS orientation)."""
        h = self.disp_shape[0]
        return dx * self.factor - 0.5, (h - dy) * self.factor - 0.5

    def _raw_to_widget(self, x, y):
        h = self.disp_shape[0]
        dx, dy = (x + 0.5) / self.factor, h - (y + 0.5) / self.factor
        o, s = self._origin(), self._scale()
        return QtCore.QPointF(o.x() + dx * s, o.y() + dy * s)

    # Painting ---------------------------------------------------------------------------
    def _stamp(self, dx, dy):
        x, y = self._disp_to_raw(dx, dy)
        gx, gy = (x + 0.5) / self.ms, (y + 0.5) / self.ms
        r = self.brush_mm / 2 / self.pixel_mm * self.factor / self.ms
        ny, nx = self.paint.shape
        y0, y1 = max(0, int(gy - r)), min(ny, int(gy + r) + 2)
        x0, x1 = max(0, int(gx - r)), min(nx, int(gx + r) + 2)
        if y0 >= y1 or x0 >= x1:
            return
        yy, xx = np.ogrid[y0:y1, x0:x1]
        disc = (xx + 0.5 - gx) ** 2 + (yy + 0.5 - gy) ** 2 <= max(r, 0.5) ** 2
        if self.tool == "paint":
            self.paint[y0:y1, x0:x1] |= disc
            self._changed.add("paint")
        elif self.tool == "protect":
            self.protect[y0:y1, x0:x1] |= disc
            self._changed.add("protect")
        elif self.tool == "texture":
            self.texture[y0:y1, x0:x1][disc] = self.texture_kind
            self._changed.add("texture")
        elif self.tool == "erase":
            self.paint[y0:y1, x0:x1] &= ~disc
            self.protect[y0:y1, x0:x1] &= ~disc
            self.texture[y0:y1, x0:x1][disc] = 0
            self._changed.update(("paint", "protect", "texture"))

    def _stroke_to(self, pt):
        d = self._to_disp(pt)
        if self._last is None:
            self._stamp(*d)
        else:
            step = max(0.5, self.brush_mm / 4 / self.pixel_mm)
            dist = np.hypot(d[0] - self._last[0], d[1] - self._last[1])
            for t in np.linspace(0, 1, max(2, int(dist / step) + 1))[1:]:
                self._stamp(self._last[0] + t * (d[0] - self._last[0]), self._last[1] + t * (d[1] - self._last[1]))
        self._last = d
        self._update_overlay()

    def mousePressEvent(self, e):
        if self.base is None or self.paint is None:
            return
        if e.button() == QtCore.Qt.MiddleButton or (e.button() == QtCore.Qt.LeftButton and self.tool == "pan"):
            self._drag = e.position()
            return
        if self.tool == "marker":
            if e.button() == QtCore.Qt.LeftButton:
                self._push_undo()
                self.markers.append(self._disp_to_raw(*self._to_disp(e.position())))
            elif e.button() == QtCore.Qt.RightButton and self.markers:
                p = e.position()
                d = [np.hypot(q.x() - p.x(), q.y() - p.y()) for q in (self._raw_to_widget(*m) for m in self.markers)]
                if min(d) < 20:
                    self._push_undo()
                    self.markers.pop(int(np.argmin(d)))
            self.update()
            self.layersChanged.emit()
            return
        if self.tool in ANNOTATION_TOOLS:
            self._annotation_press(e)
            return
        if e.button() == QtCore.Qt.LeftButton:
            self._push_undo()
            self._last = None
            self._changed = set()
            self._stroke_to(e.position())

    # Shapes and labels on the image -------------------------------------------------------
    def _raw_at(self, pt):
        x, y = self._disp_to_raw(*self._to_disp(pt))
        return [round(float(x), 2), round(float(y), 2)]

    def _annotation_press(self, e):
        if e.button() == QtCore.Qt.RightButton:
            if self._pending:
                self._finish_line()
            else:
                self._delete_annotation_at(e.position())
            return
        if e.button() != QtCore.Qt.LeftButton:
            return
        if self.tool == "label":
            hit = self._annotation_at(e.position(), kinds=("label",))
            if hit is None:
                self._edit_label(None, self._raw_at(e.position()))
        elif self.tool in ("line", "arrow"):
            self._pending.append(self._raw_at(e.position()))
            n = len(self._pending)
            self.hint.emit(f"{self.tool.capitalize()}: {n} point{'s' * (n != 1)}. Click to add more; "
                           + ("double-click or Enter to finish, Esc to cancel." if n > 1 else "Esc to cancel."))
            self.update()
        elif self.tool == "circle":
            self._circle = self._raw_at(e.position())

    def mouseDoubleClickEvent(self, e):
        if self.tool == "label" and e.button() == QtCore.Qt.LeftButton:
            hit = self._annotation_at(e.position(), kinds=("label",))
            if hit is not None:
                self._edit_label(hit, None)
        elif self.tool in ("line", "arrow") and self._pending:
            # The double-click's second press added a duplicate point.
            if len(self._pending) > 1 and self._pending[-1] == self._pending[-2]:
                self._pending.pop()
            self._finish_line()

    def keyPressEvent(self, e):
        if self._pending and e.key() in (QtCore.Qt.Key_Return, QtCore.Qt.Key_Enter):
            self._finish_line()
        elif e.key() == QtCore.Qt.Key_Escape:
            self._pending, self._circle = [], None
            self.update()
        else:
            super().keyPressEvent(e)

    def _finish_line(self):
        pts, self._pending = self._pending, []
        if len(pts) >= 2:
            self._push_undo()
            self.annotations.append({"type": "line", "points": pts, "arrow": self.tool == "arrow"})
            self.layersChanged.emit()
        self.update()

    def _edit_label(self, index, pos):
        old = self.annotations[index]["text"] if index is not None else ""
        text, ok = QtWidgets.QInputDialog.getMultiLineText(
            self, "Label on the image", "Text (one line per row). Braille cannot be shrunk, so keep it short:", old)
        if not ok:
            return
        text = "\n".join(line.strip() for line in text.split("\n") if line.strip())
        self._push_undo()
        if index is None:
            if text:
                self.annotations.append({"type": "label", "x": pos[0], "y": pos[1], "text": text})
        elif text:
            self.annotations[index]["text"] = text
        else:
            self.annotations.pop(index)
        self.layersChanged.emit()
        self.update()

    def _annotation_at(self, pt, kinds=("label", "line", "circle"), tol=12):
        """Index of the annotation under a widget point, or None."""
        best, best_d = None, tol
        s = self._scale()
        for i, a in enumerate(self.annotations):
            if a["type"] not in kinds:
                continue
            if a["type"] == "label":
                box = next((b for b, t, idx in self._label_rects() if idx == i), None)
                if box is not None and box.adjusted(-tol, -tol, tol, tol).contains(pt):
                    return i
                c = self._raw_to_widget(a["x"], a["y"])
                d = np.hypot(c.x() - pt.x(), c.y() - pt.y())
            elif a["type"] == "line":
                ws = [self._raw_to_widget(*q) for q in a["points"]]
                d = min(float(shaping._segment_distance(np.array(pt.x()), np.array(pt.y()), (p0.x(), p0.y()),
                                                        (p1.x(), p1.y()))) for p0, p1 in zip(ws[:-1], ws[1:]))
            else:
                c = self._raw_to_widget(a["x"], a["y"])
                r = a["radius_mm"] / self.pixel_mm * s
                dist = np.hypot(c.x() - pt.x(), c.y() - pt.y())
                d = max(0.0, dist - r) if a.get("filled") else abs(dist - r)
            if d < best_d:
                best, best_d = i, d
        return best

    def _delete_annotation_at(self, pt):
        i = self._annotation_at(pt)
        if i is not None:
            self._push_undo()
            self.annotations.pop(i)
            self.layersChanged.emit()
            self.update()

    def _label_rects(self):
        """(widget rect, text, annotation index) for each label, from the last preview's pads."""
        out = []
        labels = [i for i, a in enumerate(self.annotations) if a["type"] == "label"]
        if self.disp_shape is None:
            return out
        h = self.disp_shape[0]
        o, s = self._origin(), self._scale()
        for idx, i in enumerate(labels):
            a = self.annotations[i]
            if idx < len(self.label_boxes) and self.label_boxes[idx][1] == " / ".join(a["text"].split("\n")):
                c0, r0, c1, r1 = self.label_boxes[idx][0]
                rect = QtCore.QRectF(o.x() + c0 * s, o.y() + (h - r1) * s, (c1 - c0) * s, (r1 - r0) * s)
            else:   # not previewed yet
                c = self._raw_to_widget(a["x"], a["y"])
                rect = QtCore.QRectF(c.x() - 30, c.y() - 10, 60, 20)
            out.append((rect, a["text"], i))
        return out

    def mouseMoveEvent(self, e):
        self._mouse = e.position()
        if self._drag is not None:
            self.pan += e.position() - self._drag
            self._drag = e.position()
        elif self._last is not None and e.buttons() & QtCore.Qt.LeftButton:
            self._stroke_to(e.position())
        self.update()

    def _circle_radius_mm(self, pt):
        c = self._raw_to_widget(*self._circle)
        return np.hypot(pt.x() - c.x(), pt.y() - c.y()) / self._scale() * self.pixel_mm

    def mouseReleaseEvent(self, e):
        if self._drag is not None:
            self._drag = None
            return
        if self._circle is not None:
            r = self._circle_radius_mm(e.position())
            if r >= 1.0:
                self._push_undo()
                self.annotations.append({"type": "circle", "x": self._circle[0], "y": self._circle[1],
                                         "radius_mm": round(float(r), 2), "filled": self.fill_circles})
                self.layersChanged.emit()
            self._circle = None
            self.update()
            return
        if self._last is not None:
            self._last = None
            for k in self._changed:
                self.versions[k] += 1
            self._changed = set()
            self.layersChanged.emit()

    def leaveEvent(self, e):
        self._mouse = None
        self.update()

    def wheelEvent(self, e):
        if self.base is None:
            return
        p = e.position()
        before = self._to_disp(p)
        self.zoom = float(np.clip(self.zoom * (1.0015 ** e.angleDelta().y()), 0.5, 40))
        after = self._to_disp(p)
        s = self._scale()
        self.pan += QtCore.QPointF((after[0] - before[0]) * s, (after[1] - before[1]) * s)
        self.update()

    def fit(self):
        self.zoom, self.pan = 1.0, QtCore.QPointF(0, 0)
        self.update()

    def paintEvent(self, e):
        p = QtGui.QPainter(self)
        p.fillRect(self.rect(), QtGui.QColor("#1b1b1b"))
        if self.base is None:
            p.setPen(QtGui.QColor("#888"))
            p.drawText(self.rect(), QtCore.Qt.AlignCenter,
                       "Open a FITS file or image, or drag one here.\n\nHelp > Getting started explains the steps.")
            return
        h, w = self.disp_shape
        o, s = self._origin(), self._scale()
        target = QtCore.QRectF(o.x(), o.y(), w * s, h * s)
        p.setRenderHint(QtGui.QPainter.SmoothPixmapTransform, s < 2)
        p.drawImage(target, self.base)
        if self.overlay is not None:
            p.drawImage(target, self.overlay)
        p.setRenderHint(QtGui.QPainter.Antialiasing)
        r = self.marker_size_mm / 2 / self.pixel_mm * s
        for m in self.markers:
            c = self._raw_to_widget(*m)
            p.setPen(QtGui.QPen(QtGui.QColor("#ffe14d"), 2))
            p.drawEllipse(c, max(r, 4), max(r, 4))
            p.drawLine(c + QtCore.QPointF(-3, 0), c + QtCore.QPointF(3, 0))
        self._paint_annotations(p, s)
        if self._mouse is not None and self.tool in ("paint", "protect", "texture", "erase"):
            rb = self.brush_mm / 2 / self.pixel_mm * s
            p.setPen(QtGui.QPen(QtGui.QColor(255, 255, 255, 200), 1, QtCore.Qt.DashLine))
            p.drawEllipse(self._mouse, rb, rb)

    def _paint_annotations(self, p, s):
        """Shapes as they will be raised (true width), and label pads with their text."""
        width = max(1.5, self.shape_width_mm / self.pixel_mm * s)
        pen = QtGui.QPen(SHAPE_COLOUR, width, QtCore.Qt.SolidLine, QtCore.Qt.RoundCap, QtCore.Qt.RoundJoin)
        for a in self.annotations:
            if a["type"] == "line":
                pts = [self._raw_to_widget(*q) for q in a["points"]]
                p.setPen(pen)
                p.setBrush(QtCore.Qt.NoBrush)
                p.drawPolyline(QtGui.QPolygonF(pts))
                if a.get("arrow") and len(pts) > 1:
                    a0, a1 = pts[-2], pts[-1]
                    L = np.hypot(a1.x() - a0.x(), a1.y() - a0.y())
                    if L > 0:
                        ux, uy = (a1.x() - a0.x()) / L, (a1.y() - a0.y()) / L
                        hl, hw = min(3.5 * width, 0.8 * L), 1.6 * width
                        b = QtCore.QPointF(a1.x() - ux * hl, a1.y() - uy * hl)
                        head = QtGui.QPolygonF([a1, QtCore.QPointF(b.x() - uy * hw, b.y() + ux * hw),
                                                QtCore.QPointF(b.x() + uy * hw, b.y() - ux * hw)])
                        p.setPen(QtCore.Qt.NoPen)
                        p.setBrush(SHAPE_COLOUR)
                        p.drawPolygon(head)
            elif a["type"] == "circle":
                c = self._raw_to_widget(a["x"], a["y"])
                r = a["radius_mm"] / self.pixel_mm * s
                if a.get("filled"):
                    p.setPen(QtCore.Qt.NoPen)
                    p.setBrush(QtGui.QColor(255, 210, 77, 150))
                else:
                    p.setPen(pen)
                    p.setBrush(QtCore.Qt.NoBrush)
                p.drawEllipse(c, r, r)
        p.setBrush(QtCore.Qt.NoBrush)
        for rect, text, _ in self._label_rects():
            p.setPen(QtGui.QPen(LABEL_COLOUR, 1.5))
            p.fillRect(rect, QtGui.QColor(20, 40, 60, 170))
            p.drawRect(rect)
            p.setPen(QtGui.QColor("white"))
            p.drawText(rect, QtCore.Qt.AlignCenter | QtCore.Qt.TextWordWrap, text)
        # Lines and circles being drawn
        if self._pending:
            pts = [self._raw_to_widget(*q) for q in self._pending]
            if self._mouse is not None:
                pts.append(self._mouse)
            p.setPen(QtGui.QPen(SHAPE_COLOUR, width, QtCore.Qt.DashLine, QtCore.Qt.RoundCap))
            p.drawPolyline(QtGui.QPolygonF(pts))
        if self._circle is not None and self._mouse is not None:
            c = self._raw_to_widget(*self._circle)
            r = np.hypot(self._mouse.x() - c.x(), self._mouse.y() - c.y())
            p.setPen(QtGui.QPen(SHAPE_COLOUR, width, QtCore.Qt.DashLine))
            p.drawEllipse(c, r, r)
            p.setPen(QtGui.QColor("white"))
            p.drawText(self._mouse + QtCore.QPointF(10, -10), f"{r / s * self.pixel_mm:.1f} mm")


# --- Main window ----------------------------------------------------------------------------

def _spin(lo, hi, value, step=0.1, decimals=2, special=None, suffix=""):
    sb = QtWidgets.QDoubleSpinBox()
    sb.setRange(lo, hi)
    sb.setDecimals(decimals)
    sb.setSingleStep(step)
    if special:
        sb.setSpecialValueText(special)
    if suffix:
        sb.setSuffix(suffix)
    sb.setValue(value)
    return sb


def _ispin(lo, hi, value, special=None, suffix=""):
    sb = QtWidgets.QSpinBox()
    sb.setRange(lo, hi)
    if special:
        sb.setSpecialValueText(special)
    if suffix:
        sb.setSuffix(suffix)
    sb.setValue(value)
    return sb


class StretchChart(QtWidgets.QWidget):
    """What the stretch does to this image: brightness rank (faintest to brightest)
    against the height each part gets."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setMinimumHeight(160)
        self.setMinimumWidth(200)
        self.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Fixed)
        self.curve = None
        self.max_h = 10.0
        self.tactile = False
        self.edge_step = 0.2

    def set_curve(self, curve, max_height_mm, tactile, edge_step):
        self.curve, self.max_h, self.tactile, self.edge_step = curve, max_height_mm, tactile, edge_step
        self.update()

    def paintEvent(self, e):
        p = QtGui.QPainter(self)
        p.setRenderHint(QtGui.QPainter.Antialiasing)
        r = self.rect().adjusted(34, 8, -8, -30)
        p.fillRect(self.rect(), QtGui.QColor("#1d242c"))
        p.setPen(QtGui.QColor("#5c6b7a"))
        p.drawRect(r)
        small = QtGui.QFont(p.font())
        small.setPointSizeF(max(8.0, small.pointSizeF() - 2))
        p.setFont(small)
        p.setPen(QtGui.QColor("#9aa8b5"))
        p.drawText(QtCore.QRectF(r.left(), r.bottom() + 4, r.width(), 14), QtCore.Qt.AlignLeft, "faintest")
        p.drawText(QtCore.QRectF(r.left(), r.bottom() + 4, r.width(), 14), QtCore.Qt.AlignRight, "brightest")
        p.drawText(QtCore.QRectF(r.left(), r.bottom() + 15, r.width(), 14), QtCore.Qt.AlignCenter,
                   "image ranked by brightness")
        p.drawText(QtCore.QRectF(0, r.top() - 4, 32, 14), QtCore.Qt.AlignRight, f"{self.max_h:g}")
        p.drawText(QtCore.QRectF(0, r.bottom() - 10, 32, 14), QtCore.Qt.AlignRight, "0")
        p.save()
        p.translate(10, r.center().y())
        p.rotate(-90)
        p.drawText(QtCore.QRectF(-40, -8, 80, 14), QtCore.Qt.AlignCenter, "height mm")
        p.restore()
        if self.curve is None:
            p.drawText(r, QtCore.Qt.AlignCenter, "Open an image to see\nhow brightness becomes height")
            return
        if self.tactile:
            y = r.bottom() - self.edge_step * r.height()
            p.setPen(QtGui.QPen(QtGui.QColor("#ffb347"), 1, QtCore.Qt.DashLine))
            p.drawLine(QtCore.QPointF(r.left(), y), QtCore.QPointF(r.right(), y))
            p.drawText(QtCore.QRectF(r.left() + 4, y - 15, r.width(), 14), QtCore.Qt.AlignLeft, "outline step")
        flat = self.curve["flat_fraction"]
        if flat > 0.02:
            xf = r.left() + flat * r.width()
            p.fillRect(QtCore.QRectF(r.left(), r.bottom() - 6, xf - r.left(), 6), QtGui.QColor(127, 178, 229, 90))
            p.setPen(QtGui.QColor("#7fb2e5"))
            p.drawText(QtCore.QRectF(r.left() + 4, r.bottom() - 22, max(xf - r.left(), 90), 14),
                       QtCore.Qt.AlignLeft, f"{100 * flat:.0f}% flat base")
        p.setPen(QtCore.Qt.NoPen)
        p.setBrush(QtGui.QColor(255, 210, 77, 140))
        for x, h in zip(self.curve["rank"], self.curve["height"]):
            p.drawEllipse(QtCore.QPointF(r.left() + x / 100 * r.width(),
                                         r.bottom() - min(max(h, 0), 1.0) * r.height()), 1.4, 1.4)


class CollapsibleSection(QtWidgets.QWidget):
    """A settings section: a clear clickable header bar over a form that folds away."""
    HEADER_STYLE = """
        QToolButton { background: #2c3644; color: #f0f0f0; font-weight: 600; font-size: 13px;
                      border: none; border-radius: 5px; padding: 8px 10px; text-align: left; }
        QToolButton:hover { background: #374457; }
        QToolButton:checked { background: #3872a3; }
        QToolButton:checked:hover { background: #4282b8; }
    """

    def __init__(self, title, expanded=False, parent=None):
        super().__init__(parent)
        self.title = title
        self.header = QtWidgets.QToolButton()
        self.header.setText(title)
        self.header.setCheckable(True)
        self.header.setToolButtonStyle(QtCore.Qt.ToolButtonTextBesideIcon)
        self.header.setSizePolicy(QtWidgets.QSizePolicy.Ignored, QtWidgets.QSizePolicy.Fixed)
        self.header.setMinimumWidth(120)
        self.header.setStyleSheet(self.HEADER_STYLE)
        self.header.toggled.connect(self.set_expanded)
        self.help_btn = QtWidgets.QToolButton()
        self.help_btn.setText("?")
        self.help_btn.setCheckable(True)
        self.help_btn.setToolTip("What do these settings do?")
        self.help_btn.setFixedWidth(30)
        self.help_btn.setSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Preferred)
        self.help_btn.setStyleSheet(
            "QToolButton { background: #2c3644; color: #9ec9f0; font-weight: 700; border: none; border-radius: 5px; }"
            "QToolButton:hover { background: #374457; } QToolButton:checked { background: #5a7fa8; color: white; }")
        self.help_btn.toggled.connect(self._toggle_help)
        self.help_panel = QtWidgets.QLabel()
        self.help_panel.setWordWrap(True)
        self.help_panel.setTextFormat(QtCore.Qt.RichText)
        self.help_panel.setStyleSheet("background: #1f2a36; color: #dfe8f1; border-radius: 5px; padding: 8px;")
        self.help_panel.setVisible(False)
        text = help_text.SECTIONS.get(title)
        self.help_panel.setText(text or "")
        self.help_btn.setVisible(bool(text))
        self.body = QtWidgets.QFrame()
        self.body.setObjectName("sectionBody")
        self.body.setStyleSheet("#sectionBody { border-left: 2px solid #3872a3; margin-left: 6px; }")
        self.form = QtWidgets.QFormLayout(self.body)
        self.form.setContentsMargins(10, 6, 2, 10)
        self.form.setFieldGrowthPolicy(QtWidgets.QFormLayout.AllNonFixedFieldsGrow)
        self.form.setLabelAlignment(QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter)
        lay = QtWidgets.QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 4)
        lay.setSpacing(2)
        head = QtWidgets.QHBoxLayout()
        head.setSpacing(3)
        head.addWidget(self.header, 1)
        head.addWidget(self.help_btn)
        lay.addLayout(head)
        lay.addWidget(self.help_panel)
        lay.addWidget(self.body)
        self.header.setChecked(expanded)
        self.set_expanded(expanded)

    def set_expanded(self, expanded):
        if self.header.isChecked() != expanded:
            self.header.setChecked(expanded)
            return
        self.header.setArrowType(QtCore.Qt.DownArrow if expanded else QtCore.Qt.RightArrow)
        self.body.setVisible(expanded)

    def _toggle_help(self, on):
        self.help_panel.setVisible(on)

    def set_badge(self, text):
        """Short status shown on the header, e.g. 'on' or 'off'."""
        self.header.setText(f"{self.title}   ·   {text}" if text else self.title)


class AstroTouchWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("AstroTouch Desktop")
        self.resize(1400, 880)
        # ASTROTOUCH_SETTINGS_FILE keeps preferences in a file of your choice (used by the tests).
        prefs = os.environ.get("ASTROTOUCH_SETTINGS_FILE")
        self.qsettings = (QtCore.QSettings(prefs, QtCore.QSettings.IniFormat) if prefs
                          else QtCore.QSettings("AstroTouch", "AstroTouch"))

        self.fits_path = None
        self.builder = None
        self.stl_path = None
        self.last_result = None
        self.preview_worker = None
        self.preview_pending = False
        self.worker = None

        self._applying = False
        self.timer = QtCore.QTimer(self)
        self.timer.setSingleShot(True)
        self.timer.timeout.connect(self.start_preview)
        self.setAcceptDrops(True)
        self._build_menu()

        splitter = QtWidgets.QSplitter(QtCore.Qt.Horizontal)
        self.setCentralWidget(splitter)
        splitter.addWidget(self._build_sidebar())
        splitter.addWidget(self._build_views())
        splitter.setStretchFactor(1, 1)
        splitter.setSizes([460, 940])
        self._connect_changes()
        self._update_enabled()
        if self.qsettings.value("seen_intro") not in (True, "true"):
            QtCore.QTimer.singleShot(400, self.show_getting_started)
            self.qsettings.setValue("seen_intro", True)

    # Layout --------------------------------------------------------------------------
    def _build_sidebar(self):
        side = QtWidgets.QWidget()
        outer = QtWidgets.QVBoxLayout(side)
        outer.setContentsMargins(6, 6, 6, 6)

        row = QtWidgets.QHBoxLayout()
        btn_open = QtWidgets.QPushButton("Open FITS or image...")
        btn_open.clicked.connect(self.open_file)
        row.addWidget(btn_open)
        row.addWidget(QtWidgets.QLabel("HDU"))
        self.spn_hdu = _ispin(-1, 50, -1, special="auto")
        self.spn_hdu.valueChanged.connect(self.reload_file)
        row.addWidget(self.spn_hdu)
        outer.addLayout(row)
        self.lbl_file = QtWidgets.QLabel("No file selected")
        self.lbl_file.setWordWrap(True)
        outer.addWidget(self.lbl_file)

        self.sections = []
        self.sections_widget = QtWidgets.QWidget()
        self.sections_layout = QtWidgets.QVBoxLayout(self.sections_widget)
        self.sections_layout.setContentsMargins(0, 4, 6, 4)
        self.sections_layout.setSpacing(4)
        d = DEFAULTS

        # Templates
        f = self._page("Templates")
        self.cmb_plaque = QtWidgets.QComboBox()
        self.cmb_plaque.addItem("Choose a plaque format...", None)
        for key, t in pipeline.PLAQUE_FORMATS.items():
            self.cmb_plaque.addItem(t["title"], key)
            self.cmb_plaque.setItemData(self.cmb_plaque.count() - 1, t["about"], QtCore.Qt.ToolTipRole)
        self.cmb_recipe = QtWidgets.QComboBox()
        self.cmb_recipe.addItem("Choose an image recipe...", None)
        for key, t in pipeline.RECIPES.items():
            self.cmb_recipe.addItem(t["title"], key)
            self.cmb_recipe.setItemData(self.cmb_recipe.count() - 1, t["about"], QtCore.Qt.ToolTipRole)
        for cb in (self.cmb_plaque, self.cmb_recipe):
            cb.setSizeAdjustPolicy(QtWidgets.QComboBox.AdjustToMinimumContentsLengthWithIcon)
            cb.setMinimumContentsLength(14)
        self.cmb_plaque.activated.connect(lambda i: self.use_template(plaque=self.cmb_plaque.itemData(i)))
        self.cmb_recipe.activated.connect(lambda i: self.use_template(recipe=self.cmb_recipe.itemData(i)))
        self.lbl_template = self._note("")
        btns = QtWidgets.QHBoxLayout()
        b_save = QtWidgets.QPushButton("Save template...")
        b_save.clicked.connect(self.save_template)
        b_load = QtWidgets.QPushButton("Load template...")
        b_load.clicked.connect(self.load_template)
        btns.addWidget(b_save)
        btns.addWidget(b_load)
        f.addRow("Plaque format", self.cmb_plaque)
        f.addRow("Image recipe", self.cmb_recipe)
        f.addRow(self.lbl_template)
        f.addRow(btns)

        # Size and printer
        f = self._page("Model size and printer")
        self.spn_size = _spin(10, 500, d.longest_side_mm, 1, 1, suffix=" mm")
        self.spn_height = _spin(0.5, 50, d.max_height_mm, 0.5, 1, suffix=" mm")
        self.spn_base = _spin(0.5, 20, d.base_thickness_mm, 0.5, 1, suffix=" mm")
        self.spn_bed_w = _spin(50, 1000, float(self.qsettings.value("bed_w", 220)), 5, 0, suffix=" mm")
        self.spn_bed_d = _spin(50, 1000, float(self.qsettings.value("bed_d", 220)), 5, 0, suffix=" mm")
        f.addRow("Longest side of image", self.spn_size)
        f.addRow("Maximum relief", self.spn_height)
        f.addRow("Base thickness", self.spn_base)
        f.addRow("Print bed width", self.spn_bed_w)
        f.addRow("Print bed depth", self.spn_bed_d)

        # Image and stretch
        f = self._page("Image and stretch")
        self.cmb_stretch = QtWidgets.QComboBox()
        self.cmb_stretch.addItems(pipeline.STRETCHES)
        self.cmb_stretch.setCurrentText(d.stretch)
        self.cmb_stretch.setToolTip("tactile: log brightness from a definite outline step to a capped peak,\n"
                                    "with arms, knots and shells lifted. Empty sky is flat.")
        self.spn_soften = _spin(0, 1000, d.soften, 0.5, 1, special="off")
        self.spn_soften.setToolTip("Divide by N x the noise before log/asinh, so the stretch compresses\n"
                                   "whatever the data's units (needed for Jy/beam maps).")
        self.spn_clip = _spin(-1, 20, -1, 0.1, 2, special="auto", suffix=" %")
        self.spn_smooth = _spin(0, 10, d.smooth, 0.1, 2, special="off", suffix=" px")
        self.chk_invert = QtWidgets.QCheckBox("Invert heights")
        self.spn_downsample = _ispin(1, 64, d.downsample)
        self.cmb_down = QtWidgets.QComboBox()
        self.cmb_down.addItems(["mean", "linear"])
        self.spn_ridge = _spin(0, 1, d.ridge, 0.05, 2, special="off")
        self.spn_ridge.setToolTip("Mix in the block maximum so thin filaments stay continuous (cleaned maps only)")
        self.spn_min_pixel = _spin(0, 2, d.min_pixel_mm, 0.05, 2, special="off", suffix=" mm")
        f.addRow("Stretch", self.cmb_stretch)
        self.chart = StretchChart()
        f.addRow(self.chart)
        self.lbl_stretch = self._note("")
        self.lbl_stretch.setTextFormat(QtCore.Qt.RichText)
        f.addRow(self.lbl_stretch)
        f.addRow("Soften (x noise)", self.spn_soften)
        f.addRow("Clip outliers", self.spn_clip)
        f.addRow("Smoothing", self.spn_smooth)
        f.addRow("", self.chk_invert)
        f.addRow("Downsample", self.spn_downsample)
        f.addRow("Downsample method", self.cmb_down)
        f.addRow("Ridge preservation", self.spn_ridge)
        f.addRow("Smallest model pixel", self.spn_min_pixel)

        # Tactile stretch
        f = self._page("Tactile stretch")
        self.txt_edge = QtWidgets.QLineEdit()
        self.txt_edge.setPlaceholderText("auto (sky + 3 x noise)")
        self.txt_edge.setValidator(QtGui.QDoubleValidator())
        self.spn_edge_step = _spin(0, 0.9, d.edge_step, 0.05, 2)
        self.spn_contrast = _spin(0, 10, d.contrast, 0.5, 1)
        self.spn_peak_pct = _spin(90, 100, d.peak_percentile, 0.1, 2, suffix=" %")
        f.addRow("Edge (data units)", self.txt_edge)
        f.addRow("Outline step", self.spn_edge_step)
        f.addRow("Contrast boost", self.spn_contrast)
        f.addRow("Peak cap", self.spn_peak_pct)
        f.addRow(self._note("Used when the stretch is 'tactile'."))

        # Noise cleaning
        f = self._page("Noise cleaning (radio maps)")
        self.chk_clean = QtWidgets.QCheckBox("Clean noise")
        self.chk_clean.setToolTip("Mask the noise and flatten empty sky, keeping only real emission.")
        self.chk_clean.toggled.connect(self.on_clean_toggled)
        self.spn_detect = _spin(1, 20, d.detect_sigma, 0.5, 1, suffix=" σ")
        self.spn_peak_sig = _spin(1, 50, d.peak_sigma, 0.5, 1, suffix=" σ")
        self.spn_min_island = _spin(0, 1e5, 0, 1, 0, special="auto (½ beam)", suffix=" px")
        self.txt_scales = QtWidgets.QLineEdit()
        self.txt_scales.setPlaceholderText("auto (1 beam); e.g. 3 9")
        self.spn_noise_box = _ispin(0, 5000, 0, special="global noise", suffix=" px")
        self.spn_bmaj = _spin(0, 1000, 0, 0.5, 2, special="from header", suffix=" px")
        self.spn_bmin = _spin(0, 1000, 0, 0.5, 2, special="from header", suffix=" px")
        f.addRow("", self.chk_clean)
        f.addRow("Island edge", self.spn_detect)
        f.addRow("Island peak", self.spn_peak_sig)
        f.addRow("Smallest island", self.spn_min_island)
        f.addRow("Detection scales (px)", self.txt_scales)
        f.addRow("Local noise boxes", self.spn_noise_box)
        f.addRow("Beam major (FWHM)", self.spn_bmaj)
        f.addRow("Beam minor (FWHM)", self.spn_bmin)

        # Stars and sources
        f = self._page("Stars and compact sources")
        self.chk_gaia = QtWidgets.QCheckBox("Remove foreground stars (Gaia DR3)")
        self.chk_gaia.setToolTip("Needs a celestial WCS and the network. Only stars that move, or are\n"
                                 "brighter than G = 15.5, are removed; clusters inside galaxies stay.")
        self.txt_gaia_cache = QtWidgets.QLineEdit()
        self.txt_gaia_cache.setPlaceholderText("automatic (temporary folder)")
        btn_cache = QtWidgets.QPushButton("...")
        btn_cache.setFixedWidth(30)
        btn_cache.clicked.connect(self.choose_gaia_cache)
        cache_row = QtWidgets.QHBoxLayout()
        cache_row.addWidget(self.txt_gaia_cache)
        cache_row.addWidget(btn_cache)
        self.chk_points = QtWidgets.QCheckBox("Remove compact sources")
        self.chk_points.setToolTip("Faint stars, globular clusters and small background galaxies, found against\n"
                                   "a median filter. Skip for objects whose knots matter (Sc, irregulars).")
        self.spn_point_box = _ispin(3, 201, d.point_box, suffix=" px")
        self.spn_point_box.setToolTip("A few star widths (e.g. 21 for JWST pictures)")
        self.spn_point_sig = _spin(1, 50, d.point_sigma, 0.5, 1, suffix=" σ")
        self.spn_point_area = _spin(1, 1e5, d.max_point_area, 10, 0, suffix=" px")
        f.addRow("", self.chk_gaia)
        f.addRow("Gaia cache file", cache_row)
        f.addRow("", self.chk_points)
        f.addRow("Median filter width", self.spn_point_box)
        f.addRow("Detection threshold", self.spn_point_sig)
        f.addRow("Largest source", self.spn_point_area)
        f.addRow(self._note("Use the Protect brush to keep a galaxy nucleus or central star, and the "
                            "Paint-out brush for anything else (filled from its surroundings)."))

        # Paint, texture, compact islands, markers
        f = self._page("Paint-out, texture and markers")
        self.cmb_paint_mode = QtWidgets.QComboBox()
        self.cmb_paint_mode.addItems(["Fill from surroundings", "Flatten to base"])
        self.cmb_paint_mode.setToolTip("Fill: inpaint before cleaning, as if the source was never there.\n"
                                       "Flatten: set the painted area to flat base.")
        self.spn_tex_spacing = _spin(1, 20, d.texture_spacing_mm, 0.5, 1, suffix=" mm")
        self.spn_tex_height = _spin(0.1, 5, d.texture_height_mm, 0.1, 2, suffix=" mm")
        self.spn_compact = _spin(0, 50, 0, 0.5, 1, special="off", suffix=" mm")
        self.spn_compact.setToolTip("Remove isolated islands smaller than this across: they print as pips that\n"
                                    "can be mistaken for braille. Needs a flat background (cleaning or tactile).")
        self.spn_compact_frac = _spin(0, 1, d.compact_fraction, 0.05, 2, special="remove")
        self.spn_compact_cap = _spin(0, 1, d.compact_cap, 0.05, 2)
        self.spn_marker_size = _spin(0.5, 20, d.marker_size_mm, 0.5, 1, suffix=" mm")
        self.spn_marker_height = _spin(0.05, 1, d.marker_height, 0.05, 2)
        f.addRow("Paint-out", self.cmb_paint_mode)
        f.addRow("Texture spacing", self.spn_tex_spacing)
        f.addRow("Texture height", self.spn_tex_height)
        f.addRow("Remove islands under", self.spn_compact)
        f.addRow("...or keep them at", self.spn_compact_frac)
        f.addRow("...capped at", self.spn_compact_cap)
        f.addRow("Marker diameter", self.spn_marker_size)
        f.addRow("Marker height (of relief)", self.spn_marker_height)
        f.addRow(self._note("Texture dots should stay well apart from braille (2.5 mm dots); "
                            "markers are not data, so say so when presenting."))

        # Shapes and labels on the image
        f = self._page("Shapes and labels on the image")
        self.spn_shape_w = _spin(0.4, 10, d.shape_width_mm, 0.1, 1, suffix=" mm")
        self.spn_shape_h = _spin(0.1, 5, d.shape_height_mm, 0.1, 1, suffix=" mm")
        self.chk_fill_circles = QtWidgets.QCheckBox("New circles are filled discs")
        self.cmb_img_label_mode = QtWidgets.QComboBox()
        self.cmb_img_label_mode.addItems(["Braille", "Standard Text"])
        self.cmb_img_label_pad = QtWidgets.QComboBox()
        self.cmb_img_label_pad.addItems(["Flat, at base level", "Raised flat platform"])
        self.cmb_img_label_pad.setToolTip("Braille can't be read over bumpy relief, so each label sits on a flat pad:\n"
                                          "cut down to the base, or a platform at the highest point beneath it.")
        self.spn_img_label_margin = _spin(0, 10, d.image_label_margin_mm, 0.5, 1, suffix=" mm")
        f.addRow("Line and ring width", self.spn_shape_w)
        f.addRow("Line and ring height", self.spn_shape_h)
        f.addRow("", self.chk_fill_circles)
        f.addRow("Label style", self.cmb_img_label_mode)
        f.addRow("Label background", self.cmb_img_label_pad)
        f.addRow("Flat margin around label", self.spn_img_label_margin)
        f.addRow(self._note("Use the Label, Line, Arrow and Circle tools on the 2D view. Right-click a shape "
                            "or label to delete it; double-click a label to edit it. Raised lines read best "
                            "at least 1 mm wide and 0.5-1 mm high."))

        # Border and label
        f = self._page("Border and label")
        self.spn_border_w = _spin(0, 30, d.border_width_mm, 0.5, 1, suffix=" mm")
        self.spn_border_h = _spin(0, 10, d.border_height_mm, 0.5, 1, suffix=" mm")
        self.txt_label = QtWidgets.QPlainTextEdit()
        self.txt_label.setPlaceholderText("Label text, one line per row...")
        self.txt_label.setFixedHeight(self.txt_label.fontMetrics().lineSpacing() * 3 + 12)
        self.cmb_label_mode = QtWidgets.QComboBox()
        self.cmb_label_mode.addItems(["Standard Text", "Braille"])
        self.spn_label_depth = _spin(0.1, 5, DEFAULT_TEXT_DEPTH_MM, 0.1, 2, suffix=" mm")
        self.cmb_label_mode.currentTextChanged.connect(
            lambda mode: self.spn_label_depth.setValue(
                DEFAULT_BRAILLE_DOT_HEIGHT_MM if mode == "Braille" else DEFAULT_TEXT_DEPTH_MM))
        self.spn_label_size = _spin(1, 50, d.label_size_mm, 0.5, 1, suffix=" mm")
        self.spn_label_area = _spin(0, 80, d.label_area_mm, 1, 1, suffix=" mm")
        self.spn_label_spacing = _spin(0, 50, 0, 0.5, 1, special="auto", suffix=" mm")
        f.addRow("Border width", self.spn_border_w)
        f.addRow("Border height", self.spn_border_h)
        f.addRow("Label", self.txt_label)
        f.addRow("Label style", self.cmb_label_mode)
        f.addRow("Label height", self.spn_label_depth)
        f.addRow("Text cap height", self.spn_label_size)
        f.addRow("Label strip", self.spn_label_area)
        f.addRow("Line spacing", self.spn_label_spacing)

        self.sections_layout.addStretch(1)
        fold_row = QtWidgets.QHBoxLayout()
        fold_row.addWidget(QtWidgets.QLabel("<b>Settings</b>"))
        fold_row.addStretch(1)
        for text, state in (("Expand all", True), ("Collapse all", False)):
            b = QtWidgets.QPushButton(text)
            b.setFlat(True)
            b.setStyleSheet("color: #7fb2e5;")
            b.clicked.connect(lambda _=False, st=state: [sec.set_expanded(st) for sec in self.sections])
            fold_row.addWidget(b)
        outer.addLayout(fold_row)
        scroll = QtWidgets.QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QtWidgets.QFrame.NoFrame)
        scroll.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        scroll.setWidget(self.sections_widget)
        outer.addWidget(scroll, 1)

        self.lbl_summary = QtWidgets.QLabel("")
        self.lbl_summary.setWordWrap(True)
        self.lbl_summary.setTextFormat(QtCore.Qt.RichText)
        outer.addWidget(self.lbl_summary)

        self.btn_process = QtWidgets.QPushButton("PROCESS 3D MODEL")
        self.btn_process.setFixedHeight(46)
        self.btn_process.setStyleSheet("background-color: #3872a3; color: white; font-weight: bold;")
        self.btn_process.clicked.connect(self.start_processing)
        outer.addWidget(self.btn_process)
        self.progress_bar = QtWidgets.QProgressBar()
        self.progress_bar.setVisible(False)
        self.progress_bar.setRange(0, 0)
        outer.addWidget(self.progress_bar)
        self.btn_save = QtWidgets.QPushButton("SAVE STL AS...")
        self.btn_save.clicked.connect(self.save_stl)
        self.btn_save.setEnabled(False)
        outer.addWidget(self.btn_save)
        self.status = QtWidgets.QLabel("Ready")
        self.status.setWordWrap(True)
        outer.addWidget(self.status)
        side.setMinimumWidth(380)
        return side

    OPEN_BY_DEFAULT = ("Templates", "Image and stretch")

    def _page(self, title):
        saved = self.qsettings.value(f"sections/{title}")
        expanded = (saved in (True, "true")) if saved is not None else title in self.OPEN_BY_DEFAULT
        sec = CollapsibleSection(title, expanded)
        sec.header.toggled.connect(lambda on, t=title: self.qsettings.setValue(f"sections/{t}", on))
        self.sections.append(sec)
        self.sections_layout.addWidget(sec)
        return sec.form

    def _section(self, title):
        return next(sec for sec in self.sections if sec.title == title)

    def _update_badges(self, *args):
        chosen = [table[cb.currentData()]["short"] for cb, table in ((self.cmb_plaque, pipeline.PLAQUE_FORMATS),
                                                                     (self.cmb_recipe, pipeline.RECIPES))
                  if cb.currentData() is not None]
        self._section("Templates").set_badge(" + ".join(chosen) if chosen else "none")
        stretch = self.cmb_stretch.currentText()
        self._section("Image and stretch").set_badge(stretch)
        self._section("Tactile stretch").set_badge("in use" if stretch == "tactile" else "off")
        self._section("Noise cleaning (radio maps)").set_badge("on" if self.chk_clean.isChecked() else "off")
        stars = [n for n, on in (("Gaia", self.chk_gaia.isChecked()), ("compact", self.chk_points.isChecked())) if on]
        self._section("Stars and compact sources").set_badge(" + ".join(stars) if stars else "off")
        anns = self.canvas.annotations
        n_lab = sum(a["type"] == "label" for a in anns)
        n_shape = len(anns) - n_lab
        self._section("Shapes and labels on the image").set_badge(
            ", ".join(t for t in (f"{n_shape} shape{'s' * (n_shape != 1)}" if n_shape else "",
                                  f"{n_lab} label{'s' * (n_lab != 1)}" if n_lab else "") if t) or "none")
        labels = [line for line in self.txt_label.toPlainText().split("\n") if line.strip()]
        self._section("Border and label").set_badge(
            f"{len(labels)} line{'s' * (len(labels) != 1)}, {self.cmb_label_mode.currentText().lower()}"
            if labels else "no label")

    @staticmethod
    def _note(text):
        lbl = QtWidgets.QLabel(text)
        lbl.setWordWrap(True)
        lbl.setStyleSheet("color: #999; font-size: 11px;")
        return lbl

    def _build_views(self):
        self.tabs = QtWidgets.QTabWidget()

        view2d = QtWidgets.QWidget()
        lay = QtWidgets.QVBoxLayout(view2d)
        lay.setContentsMargins(0, 0, 0, 0)
        bar = QtWidgets.QToolBar()
        self.tool_group = QtGui.QActionGroup(self)
        self.tool_actions = {}
        for key, name, tip in TOOLS:
            a = QtGui.QAction(name, self)
            a.setCheckable(True)
            a.setToolTip(tip)
            a.triggered.connect(lambda _=False, k=key: self.set_tool(k))
            self.tool_group.addAction(a)
            bar.addAction(a)
            self.tool_actions[key] = a
        bar2 = QtWidgets.QToolBar()
        bar2.addWidget(QtWidgets.QLabel(" Brush "))
        self.spn_brush = _spin(0.2, 100, 3.0, 0.5, 1, suffix=" mm")
        self.spn_brush.valueChanged.connect(self._brush_changed)
        bar2.addWidget(self.spn_brush)
        bar2.addWidget(QtWidgets.QLabel(" Texture "))
        self.cmb_texture = QtWidgets.QComboBox()
        self.cmb_texture.addItems([t.capitalize() for t in shaping.TEXTURES])
        self.cmb_texture.currentIndexChanged.connect(self._texture_kind_changed)
        bar2.addWidget(self.cmb_texture)
        bar2.addSeparator()
        self.cmb_display = QtWidgets.QComboBox()
        self.cmb_display.addItems(["Shaded relief", "Input image"])
        self.cmb_display.currentIndexChanged.connect(self.show_result)
        bar2.addWidget(self.cmb_display)
        self.chk_show_removed = QtWidgets.QCheckBox("Show removed")
        self.chk_show_removed.setChecked(True)
        self.chk_show_removed.setToolTip("Outline what star and source removal took out:\n"
                                         "cyan = Gaia stars, magenta = compact sources.\n"
                                         "Switch to 'Input image' to see the originals.")
        self.chk_show_removed.toggled.connect(lambda on: self.canvas.set_show_removed(on))
        bar2.addWidget(self.chk_show_removed)
        a_undo = QtGui.QAction("Undo", self)
        a_undo.setShortcut(QtGui.QKeySequence.Undo)
        a_undo.triggered.connect(lambda: self.canvas.undo())
        bar2.addAction(a_undo)
        clear_btn = QtWidgets.QToolButton()
        clear_btn.setText("Clear")
        clear_btn.setPopupMode(QtWidgets.QToolButton.InstantPopup)
        menu = QtWidgets.QMenu(clear_btn)
        for label, what in (("Paint-out", "paint"), ("Protect", "protect"), ("Texture", "texture"),
                            ("Markers", "markers"), ("Shapes and labels", "annotations")):
            menu.addAction(label, lambda w=what: self.canvas.clear(w))
        clear_btn.setMenu(menu)
        bar2.addWidget(clear_btn)
        a_fit = QtGui.QAction("Fit", self)
        a_fit.triggered.connect(lambda: self.canvas.fit())
        bar2.addAction(a_fit)
        lay.addWidget(bar)
        lay.addWidget(bar2)
        self.lbl_hint = QtWidgets.QLabel(help_text.TOOLS["paint"])
        self.lbl_hint.setStyleSheet("background: #22303d; color: #cfe3f5; padding: 5px 10px;")
        self.lbl_hint.setWordWrap(True)
        lay.addWidget(self.lbl_hint)
        self.canvas = PaintCanvas()
        self.canvas.hint.connect(self.lbl_hint.setText)
        self.canvas.layersChanged.connect(self.schedule_preview)
        lay.addWidget(self.canvas, 1)
        self.tabs.addTab(view2d, "2D view and painting")

        self.plotter_widget = QtWidgets.QWidget()
        pl = QtWidgets.QVBoxLayout(self.plotter_widget)
        pl.setContentsMargins(0, 0, 0, 0)
        self.plotter = QtInteractor(self.plotter_widget)
        pl.addWidget(self.plotter)
        self.plotter.set_background("black")
        self.plotter.enable_eye_dome_lighting()
        btn_recenter = QtWidgets.QPushButton("Center camera")
        btn_recenter.clicked.connect(self.recenter_view)
        pl.addWidget(btn_recenter)
        self.tabs.addTab(self.plotter_widget, "3D model")

        self.tool_actions["paint"].setChecked(True)
        return self.tabs

    def _connect_changes(self):
        """Any setting change rebuilds the preview a moment later."""
        for w in self.findChildren(QtWidgets.QAbstractSpinBox):
            if w is self.spn_hdu or w is self.spn_brush:
                continue
            w.valueChanged.connect(self.schedule_preview)
        for w in (self.cmb_stretch, self.cmb_down, self.cmb_paint_mode, self.cmb_label_mode):
            w.currentIndexChanged.connect(self.schedule_preview)
        for w in (self.chk_invert, self.chk_clean, self.chk_gaia, self.chk_points):
            w.toggled.connect(self.schedule_preview)
        for w in (self.txt_edge, self.txt_scales, self.txt_gaia_cache):
            w.editingFinished.connect(self.schedule_preview)
        self.txt_label.textChanged.connect(self.schedule_preview)
        self.cmb_stretch.currentIndexChanged.connect(self._update_enabled)
        for w in (self.cmb_stretch, self.cmb_label_mode):
            w.currentIndexChanged.connect(self._update_badges)
        for w in (self.chk_clean, self.chk_gaia, self.chk_points):
            w.toggled.connect(self._update_badges)
        self.txt_label.textChanged.connect(self._update_badges)
        self._update_badges()
        self.spn_marker_size.valueChanged.connect(self._marker_size_changed)
        self.spn_shape_w.valueChanged.connect(lambda v: (setattr(self.canvas, "shape_width_mm", v),
                                                         self.canvas.update()))
        self.chk_fill_circles.toggled.connect(lambda on: setattr(self.canvas, "fill_circles", on))
        for w in (self.cmb_img_label_mode, self.cmb_img_label_pad):
            w.currentIndexChanged.connect(self.schedule_preview)
        self.canvas.layersChanged.connect(self._update_badges)
        self.spn_bed_w.valueChanged.connect(lambda v: self.qsettings.setValue("bed_w", v))
        self.spn_bed_d.valueChanged.connect(lambda v: self.qsettings.setValue("bed_d", v))

    def _update_enabled(self):
        tactile = self.cmb_stretch.currentText() == "tactile"
        for w in (self.txt_edge, self.spn_edge_step, self.spn_contrast, self.spn_peak_pct):
            w.setEnabled(tactile)
        self.spn_soften.setEnabled(self.cmb_stretch.currentText() in ("log", "asinh"))

    # Tools ----------------------------------------------------------------------------------
    def set_tool(self, key):
        self.lbl_hint.setText(help_text.TOOLS.get(key, ""))
        self.canvas.tool = key
        self.tool_actions[key].setChecked(True)
        self.canvas._pending, self.canvas._circle = [], None
        self.canvas.setCursor(QtCore.Qt.OpenHandCursor if key == "pan" else
                              QtCore.Qt.IBeamCursor if key == "label" else
                              QtCore.Qt.CrossCursor if key in ("marker", "line", "arrow", "circle")
                              else QtCore.Qt.BlankCursor)
        self.canvas.setFocus()
        self.canvas.update()

    def _brush_changed(self, v):
        self.canvas.brush_mm = v

    def _texture_kind_changed(self, i):
        self.canvas.texture_kind = i + 1
        self.set_tool("texture")

    def _marker_size_changed(self, v):
        self.canvas.marker_size_mm = v
        self.canvas.update()

    # Settings ---------------------------------------------------------------------------------
    def settings(self):
        def opt(sb):
            return sb.value() if sb.value() != sb.minimum() else None

        edge = self.txt_edge.text().strip()
        scales = tuple(float(v) for v in self.txt_scales.text().replace(",", " ").split()) or None
        beam = (self.spn_bmaj.value(), self.spn_bmin.value() or self.spn_bmaj.value()) \
            if self.spn_bmaj.value() > 0 else None
        cache = self.txt_gaia_cache.text().strip() or None
        if self.chk_gaia.isChecked() and cache is None and self.fits_path:
            tag = hashlib.sha1(f"{self.fits_path}:{self.spn_hdu.value()}".encode()).hexdigest()[:12]
            cache = str(Path(tempfile.gettempdir()) / f"astrotouch_gaia_{tag}.csv")
        braille = self.cmb_label_mode.currentText() == "Braille"
        return pipeline.Settings(
            longest_side_mm=self.spn_size.value(), max_height_mm=self.spn_height.value(),
            base_thickness_mm=self.spn_base.value(), clip=opt(self.spn_clip), smooth=self.spn_smooth.value(),
            downsample=self.spn_downsample.value(), downsample_method=self.cmb_down.currentText(),
            ridge=self.spn_ridge.value(), min_pixel_mm=self.spn_min_pixel.value(),
            stretch=self.cmb_stretch.currentText(), soften=self.spn_soften.value(),
            invert=self.chk_invert.isChecked(), edge=float(edge) if edge else None,
            edge_step=self.spn_edge_step.value(), contrast=self.spn_contrast.value(),
            peak_percentile=self.spn_peak_pct.value(), clean=self.chk_clean.isChecked(), beam=beam,
            detect_sigma=self.spn_detect.value(), peak_sigma=self.spn_peak_sig.value(),
            min_island=opt(self.spn_min_island), smooth_scales=scales,
            noise_box=self.spn_noise_box.value() or None, mask_stars=self.chk_gaia.isChecked(),
            gaia_cache=cache, remove_points=self.chk_points.isChecked(), point_box=self.spn_point_box.value(),
            point_sigma=self.spn_point_sig.value(), max_point_area=self.spn_point_area.value(),
            paint_mode="fill" if self.cmb_paint_mode.currentIndex() == 0 else "flatten",
            remove_compact=opt(self.spn_compact), compact_fraction=self.spn_compact_frac.value(),
            compact_cap=self.spn_compact_cap.value(), markers=tuple(self.canvas.markers),
            marker_size_mm=self.spn_marker_size.value(), marker_height=self.spn_marker_height.value(),
            texture_spacing_mm=self.spn_tex_spacing.value(), texture_height_mm=self.spn_tex_height.value(),
            annotations=tuple(dict(a) for a in self.canvas.annotations),
            shape_width_mm=self.spn_shape_w.value(), shape_height_mm=self.spn_shape_h.value(),
            image_label_mode="braille" if self.cmb_img_label_mode.currentIndex() == 0 else "text",
            image_label_pad="base" if self.cmb_img_label_pad.currentIndex() == 0 else "raised",
            image_label_margin_mm=self.spn_img_label_margin.value(),
            border_width_mm=self.spn_border_w.value(), border_height_mm=self.spn_border_h.value(),
            label=tuple(line for line in self.txt_label.toPlainText().split("\n") if line.strip()),
            label_mode="braille" if braille else "text", label_depth_mm=self.spn_label_depth.value(),
            label_size_mm=self.spn_label_size.value(), label_area_mm=self.spn_label_area.value(),
            label_spacing_mm=opt(self.spn_label_spacing))

    def on_clean_toggled(self, checked):
        if checked and self.cmb_stretch.currentText() == "linear":
            self.cmb_stretch.setCurrentText("asinh")

    def choose_gaia_cache(self):
        path, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Gaia cache file", "", "CSV (*.csv)")
        if path:
            self.txt_gaia_cache.setText(path)
            self.schedule_preview()

    # Files --------------------------------------------------------------------------------------
    def open_file(self):
        path, _ = QtWidgets.QFileDialog.getOpenFileName(
            self, "Open FITS or image", "",
            "FITS or image (*.fits *.fit *.fz *.png *.jpg *.jpeg *.tif *.tiff);;FITS (*.fits *.fit *.fz);;"
            "Images (*.png *.jpg *.jpeg *.tif *.tiff)")
        if path:
            self.load_file(path)

    def load_file(self, path):
        self.fits_path = path
        self.lbl_file.setText(os.path.basename(path))
        self.reload_file()

    def reload_file(self):
        if not self.fits_path:
            return
        hdu = self.spn_hdu.value()
        try:
            self.builder = pipeline.ModelBuilder.from_file(self.fits_path, None if hdu < 0 else hdu)
        except Exception as e:
            self.builder = None
            self.set_status(f"Could not open: {e}", error=True)
            return
        if self.canvas.raw_shape != self.builder.raw.shape:
            self.canvas.reset_layers(self.builder.raw.shape)
        ny, nx = self.builder.raw.shape
        self.lbl_file.setText(f"{os.path.basename(self.fits_path)} ({nx} x {ny} px)")
        self.tabs.setCurrentIndex(0)
        self.schedule_preview()

    # Preview ----------------------------------------------------------------------------------
    def schedule_preview(self, *args):
        if self.builder is not None and not self._applying:
            self.timer.start(PREVIEW_DELAY_MS)

    def start_preview(self):
        if self.builder is None:
            return
        if self.preview_worker is not None and self.preview_worker.isRunning():
            self.preview_pending = True
            return
        self.preview_pending = False
        self.set_status("Updating preview...")
        self.preview_worker = PreviewWorker(self.builder, self.settings(), self.canvas.layers())
        self.preview_worker.done.connect(self.on_preview)
        self.preview_worker.failed.connect(self.on_preview_failed)
        self.preview_worker.finished.connect(self._preview_finished)
        self.preview_worker.start()

    def _preview_finished(self):
        if self.preview_pending:
            self.start_preview()

    def on_preview(self, result, settings):
        self.last_result = (result, settings)
        self.show_result()
        size = pipeline.plaque_size(settings, result)
        bed = (self.spn_bed_w.value(), self.spn_bed_d.value())
        fits_ok = pipeline.fits_bed(size, bed)
        bed_txt = (f"fits the {bed[0]:.0f} x {bed[1]:.0f} mm bed" if fits_ok else
                   f"<b style='color:#ff6b6b'>does not fit the {bed[0]:.0f} x {bed[1]:.0f} mm bed</b>")
        self.lbl_summary.setText(
            f"Plaque ≈ <b>{size[0]:.0f} x {size[1]:.0f} x {size[2]:.1f} mm</b>, {bed_txt}.<br>"
            f"{result.heights.shape[1]} x {result.heights.shape[0]} model pixels of {result.pixel_mm:.2f} mm"
            + (f" (downsampled x{result.factor})" if result.factor > 1 else "")
            + self._removed_summary(settings, result))
        if result.warnings:
            self.set_status("Check: " + " ".join(result.warnings), warning=True)
        else:
            self.set_status("Preview up to date.")

    @staticmethod
    def _removed_summary(settings, result):
        r = result.removed
        if r is None:
            return ""
        parts = []
        if settings.mask_stars:
            parts.append(f"<span style='color:#00e6ff'>{r['n_stars']} Gaia star{'s' * (r['n_stars'] != 1)}</span>")
        if settings.remove_points:
            parts.append(f"<span style='color:#ff46dc'>{r['n_points']} compact "
                         f"source{'s' * (r['n_points'] != 1)}</span>")
        if not parts:
            return ""
        return f"<br>Removed {' and '.join(parts)} ({100 * r['fraction']:.1f}% of the image)." 

    def on_preview_failed(self, message):
        self.set_status(f"Preview: {message}", error=True)

    def show_result(self, *args):
        if self.last_result is None:
            return
        result, settings = self.last_result
        if self.cmb_display.currentIndex() == 1 and result.display is not None:
            img = np.flipud(result.display)
        else:
            h = np.flipud(np.clip(result.heights, 0, None))
            img = 0.7 * shaping.hillshade(h, result.pixel_mm, settings.max_height_mm) + 0.3 * np.clip(h, 0, 1)
        self._update_chart(result, settings)
        self.canvas.label_boxes = [(box, text) for _, _, _, _, box, text in result.image_labels]
        self.canvas.set_image(img, result.factor, result.pixel_mm)
        self.canvas.set_removed(result.removed)

    def set_status(self, text, error=False, warning=False):
        self.status.setText(text)
        self.status.setStyleSheet("color: #ff6b6b;" if error else "color: #ffb347;" if warning else "")

    # Build and save ----------------------------------------------------------------------------
    def start_processing(self):
        if self.builder is None:
            return
        self.btn_process.setEnabled(False)
        self.progress_bar.setVisible(True)
        self.worker = BuildWorker(self.builder, self.settings(), self.canvas.layers())
        self.worker.progress.connect(self.set_status)
        self.worker.done.connect(self.update_view)
        self.worker.failed.connect(self.handle_error)
        self.worker.start()

    def handle_error(self, message):
        QtWidgets.QMessageBox.critical(self, "Error", message)
        self.btn_process.setEnabled(True)
        self.progress_bar.setVisible(False)
        self.set_status(message, error=True)

    def update_view(self, result):
        self.stl_path = result['path']
        self.btn_process.setEnabled(True)
        self.btn_save.setEnabled(True)
        self.progress_bar.setVisible(False)
        w, d, h = result['size']
        self.set_status(f"Done: {w:.0f} x {d:.0f} x {h:.1f} mm, {result['faces']:,} faces.")
        self.plotter.clear()
        mesh = pv.read(str(self.stl_path))
        # Centre the mesh for the view only; the exported STL is unchanged.
        mesh.translate(-np.array(mesh.center), inplace=True)
        self.plotter.add_mesh(mesh, color="lightgray", smooth_shading=True, name="model")
        self.tabs.setCurrentIndex(1)
        self.recenter_view(mesh=mesh)

    def recenter_view(self, *args, mesh=None):
        """Look straight down at the model."""
        if mesh is None:
            try:
                for actor in self.plotter.renderer.actors.values():
                    if getattr(actor, "mapper", None) is not None:
                        mesh = actor.mapper.dataset
                        break
            except Exception:
                pass

        def do_recenter():
            self.plotter.view_xy()
            if mesh is not None and hasattr(mesh, 'bounds'):
                self.plotter.reset_camera(bounds=mesh.bounds)
            else:
                self.plotter.reset_camera()
            self.plotter.render()
        QtCore.QTimer.singleShot(100, do_recenter)

    def save_stl(self):
        if not self.stl_path:
            return
        save_path, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Save STL", "", "STL (*.stl)")
        if save_path:
            saved = self.save_outputs(save_path)
            QtWidgets.QMessageBox.information(self, "Saved", "Saved:\n" + "\n".join(saved))

    def save_outputs(self, save_path):
        """Saves the STL, plus the painted layers and the settings beside it, so the
        model can be rebuilt (also from the command line)."""
        save_path = Path(save_path)
        shutil.copy(self.stl_path, save_path)
        saved = [save_path.name]
        stem = save_path.with_suffix("")
        layers = {}
        for name, arr, flag in (("paintout", self.canvas.paint, "--paint-out"),
                                ("protect", self.canvas.protect, "--protect-mask"),
                                ("texture", self.canvas.texture, "--texture-map")):
            if arr is not None and arr.any():
                path = Path(f"{stem}.{name}.fits")
                fits.writeto(path, arr.astype(np.uint8), overwrite=True)
                layers[flag] = path.name
                saved.append(path.name)
        meta = {"input": self.fits_path, "hdu": None if self.spn_hdu.value() < 0 else self.spn_hdu.value(),
                "settings": self.settings().to_dict(), "layers": layers,
                "note": "Painted layers are on a grid covering the whole input image; pass them to "
                        "cli.py with the flags given as keys above."}
        if self.canvas.annotations:
            layers["--annotations"] = f"{stem.name}.astrotouch.json"
        meta_path = Path(f"{stem}.astrotouch.json")
        meta_path.write_text(json.dumps(meta, indent=2))
        saved.append(meta_path.name)
        return saved

    # Help ----------------------------------------------------------------------------------------
    def _build_menu(self):
        m = self.menuBar().addMenu("&Help")
        m.addAction("Getting started", self.show_getting_started)
        a = m.addAction("Show tool hints")
        a.setCheckable(True)
        a.setChecked(True)
        a.toggled.connect(lambda on: self.lbl_hint.setVisible(on))
        m.addAction("Explain every section", lambda: [sec.help_btn.setChecked(True) or sec.set_expanded(True)
                                                       for sec in self.sections if sec.help_btn.isVisible()])
        m.addAction("Hide all explanations", lambda: [sec.help_btn.setChecked(False) for sec in self.sections])

    def show_getting_started(self):
        dlg = QtWidgets.QDialog(self)
        dlg.setWindowTitle("Getting started with AstroTouch")
        dlg.resize(620, 620)
        lay = QtWidgets.QVBoxLayout(dlg)
        view = QtWidgets.QTextBrowser()
        view.setHtml(help_text.GETTING_STARTED)
        view.setOpenExternalLinks(True)
        lay.addWidget(view)
        ok = QtWidgets.QPushButton("Close")
        ok.clicked.connect(dlg.accept)
        lay.addWidget(ok, alignment=QtCore.Qt.AlignRight)
        dlg.setModal(False)
        dlg.show()
        self._intro = dlg

    def dragEnterEvent(self, e):
        urls = e.mimeData().urls() if e.mimeData().hasUrls() else []
        if urls and urls[0].toLocalFile().lower().endswith(OPENABLE):
            e.acceptProposedAction()

    def dropEvent(self, e):
        self.load_file(e.mimeData().urls()[0].toLocalFile())

    # Templates ---------------------------------------------------------------------------------
    def use_template(self, plaque=None, recipe=None, custom=None):
        if plaque is None and recipe is None and custom is None:
            return
        self.apply_settings(pipeline.apply_template(self.settings(), plaque=plaque, recipe=recipe, custom=custom))
        parts = []
        if plaque:
            parts.append(pipeline.PLAQUE_FORMATS[plaque]["about"])
        if recipe:
            parts.append(pipeline.RECIPES[recipe]["about"])
        if parts:
            self.lbl_template.setText(" ".join(parts))
        self._update_badges()

    def save_template(self):
        path, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Save template", "my_template.astrotouch-template.json",
                                                        "Template (*.json)")
        if path:
            Path(path).write_text(json.dumps({"astrotouch_template": 1,
                                              "settings": pipeline.template_dict(self.settings())}, indent=2))
            self.set_status(f"Template saved: {os.path.basename(path)}")

    def load_template(self):
        path, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Load template", "", "Template (*.json)")
        if not path:
            return
        try:
            data = json.loads(Path(path).read_text())
            custom = data.get("settings", data)
            self.use_template(custom={k: v for k, v in custom.items() if k not in pipeline.TEMPLATE_EXCLUDED})
            self.lbl_template.setText(f"Loaded {os.path.basename(path)}.")
        except Exception as e:
            self.set_status(f"Could not load template: {e}", error=True)

    def apply_settings(self, s):
        """Shows a Settings object in the widgets (everything except per-image content)."""
        self._applying = True
        try:
            def setv(sb, v, none=None):
                sb.setValue(sb.minimum() if v is None and none is None else (none if v is None else v))
            setv(self.spn_size, s.longest_side_mm)
            setv(self.spn_height, s.max_height_mm)
            setv(self.spn_base, s.base_thickness_mm)
            self.chk_clean.setChecked(bool(s.clean))
            self.cmb_stretch.setCurrentText(s.stretch)
            setv(self.spn_soften, s.soften)
            setv(self.spn_clip, s.clip)
            setv(self.spn_smooth, s.smooth)
            self.chk_invert.setChecked(bool(s.invert))
            setv(self.spn_downsample, s.downsample)
            self.cmb_down.setCurrentText(s.downsample_method)
            setv(self.spn_ridge, s.ridge)
            setv(self.spn_min_pixel, s.min_pixel_mm)
            self.txt_edge.setText("" if s.edge is None else f"{s.edge:g}")
            setv(self.spn_edge_step, s.edge_step)
            setv(self.spn_contrast, s.contrast)
            setv(self.spn_peak_pct, s.peak_percentile)
            setv(self.spn_bmaj, s.beam[0] if s.beam else 0)
            setv(self.spn_bmin, s.beam[1] if s.beam else 0)
            setv(self.spn_detect, s.detect_sigma)
            setv(self.spn_peak_sig, s.peak_sigma)
            setv(self.spn_min_island, s.min_island, 0)
            self.txt_scales.setText(" ".join(f"{v:g}" for v in s.smooth_scales) if s.smooth_scales else "")
            setv(self.spn_noise_box, s.noise_box, 0)
            self.chk_gaia.setChecked(bool(s.mask_stars))
            self.chk_points.setChecked(bool(s.remove_points))
            setv(self.spn_point_box, s.point_box)
            setv(self.spn_point_sig, s.point_sigma)
            setv(self.spn_point_area, s.max_point_area)
            self.cmb_paint_mode.setCurrentIndex(0 if s.paint_mode == "fill" else 1)
            setv(self.spn_compact, s.remove_compact, 0)
            setv(self.spn_compact_frac, s.compact_fraction)
            setv(self.spn_compact_cap, s.compact_cap)
            setv(self.spn_marker_size, s.marker_size_mm)
            setv(self.spn_marker_height, s.marker_height)
            setv(self.spn_tex_spacing, s.texture_spacing_mm)
            setv(self.spn_tex_height, s.texture_height_mm)
            setv(self.spn_shape_w, s.shape_width_mm)
            setv(self.spn_shape_h, s.shape_height_mm)
            self.cmb_img_label_mode.setCurrentIndex(0 if s.image_label_mode == "braille" else 1)
            self.cmb_img_label_pad.setCurrentIndex(0 if s.image_label_pad == "base" else 1)
            setv(self.spn_img_label_margin, s.image_label_margin_mm)
            setv(self.spn_border_w, s.border_width_mm)
            setv(self.spn_border_h, s.border_height_mm)
            self.cmb_label_mode.setCurrentText("Braille" if s.label_mode == "braille" else "Standard Text")
            if s.label_depth_mm is not None:
                setv(self.spn_label_depth, s.label_depth_mm)
            setv(self.spn_label_size, s.label_size_mm)
            setv(self.spn_label_area, s.label_area_mm)
            setv(self.spn_label_spacing, s.label_spacing_mm, 0)
        finally:
            self._applying = False
        self._update_enabled()
        self._update_badges()
        self.schedule_preview()

    # Stretch chart -------------------------------------------------------------------------------
    STRETCH_NOTES = {
        "linear": "Height is in proportion to brightness.",
        "log": "Logarithmic: bright parts are compressed so fainter ones get height too.",
        "asinh": "asinh: like log for bright parts, but gentle and linear near zero, so noise is not exaggerated.",
        "tactile": "Tactile: flat below the edge, a step at the outline, log brightness up to a capped summit; "
                   "points above the curve were lifted by the contrast boost.",
    }

    def _update_chart(self, result, settings):
        c = result.curve
        self.chart.set_curve(c, settings.max_height_mm, settings.stretch == "tactile", settings.edge_step)
        if c is None:
            self.lbl_stretch.setText("")
            return
        flat, top = 100 * c["flat_fraction"], 100 * c["top_fraction"]
        msg = [self.STRETCH_NOTES.get(settings.stretch, "")]
        msg.append(f"<b>{flat:.0f}%</b> of the image is flat base; <b>{top:.1f}%</b> reaches the top.")
        if flat < 20 and settings.stretch in ("log", "asinh", "linear") and not settings.clean:
            msg.append("Little is flat, so noise will be felt as texture: try noise cleaning or the tactile stretch.")
        if top > 15:
            msg.append("Much of the image is at full height: lower the stretch or raise the peak cap.")
        self.lbl_stretch.setText(" ".join(msg))


def main():
    app = QtWidgets.QApplication(sys.argv)
    app.setStyle("Fusion")
    app.setApplicationName("AstroTouch")
    window = AstroTouchWindow()
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
