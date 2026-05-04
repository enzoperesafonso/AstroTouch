#!/usr/bin/env python3
import sys
import os
from pathlib import Path
import tempfile

import numpy as np
from astropy.io import fits
from PySide6 import QtWidgets, QtCore, QtGui
import pyvista as pv
from pyvistaqt import QtInteractor

from scipy.ndimage import gaussian_filter

# Import core logic
try:
    from .core import load_fits_data, preprocess_image, generate_mesh
except ImportError:
    from core import load_fits_data, preprocess_image, generate_mesh

class Worker(QtCore.QThread):
    finished = QtCore.Signal(object)
    error = QtCore.Signal(str)
    progress = QtCore.Signal(str)

    def __init__(self, params):
        super().__init__()
        self.params = params

    def run(self):
        try:
            self.progress.emit("Loading FITS data...")
            data = load_fits_data(self.params['fits_path'], self.params['hdu'])
            
            self.progress.emit("Preprocessing image...")
            norm_data = preprocess_image(
                data,
                downsample_factor=self.params['downsample'],
                clip_percentile=self.params['clip'],
                log_scale=(self.params['scale_mode'] == 'log'),
                asinh_scale=(self.params['scale_mode'] == 'asinh'),
                invert=self.params['invert']
            )
            
            self.progress.emit("Generating 3D mesh...")
            stl_mesh = generate_mesh(
                norm_data,
                longest_side_mm=self.params['longest_side'],
                max_height_mm=self.params['max_height'],
                base_thickness_mm=self.params['base_thickness'],
                smoothing_sigma=self.params['smooth'],
                border_width_mm=self.params['border_width'],
                border_height_mm=self.params['border_height'],
                label_text=self.params.get('label_text', ""),
                label_mode=self.params.get('label_mode', "text"),
                label_depth_mm=self.params.get('label_depth', 1.5),
                label_font_size_mm=self.params.get('label_size', 5.0),
                label_area_mm=self.params.get('label_area', 0.0)
            )
            
            temp_stl = Path(tempfile.gettempdir()) / "preview.stl"
            stl_mesh.save(str(temp_stl))
            
            self.finished.emit({
                'path': temp_stl,
                'vertices': len(stl_mesh.vectors) * 3,
                'faces': len(stl_mesh.vectors)
            })
        except Exception as e:
            self.error.emit(str(e))

class AstroTouchWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("AstroTouch Desktop")
        self.resize(1280, 800)

        self.fits_path = None
        self.stl_path = None
        self.raw_data = None
        self.last_hdu = -1

        # --- Layout using QSplitter for proper expansion ---
        splitter = QtWidgets.QSplitter(QtCore.Qt.Horizontal)
        self.setCentralWidget(splitter)

        # Sidebar Container with Scroll Area
        self.scroll_area = QtWidgets.QScrollArea()
        self.scroll_area.setWidgetResizable(True)
        self.scroll_area.setFrameShape(QtWidgets.QFrame.NoFrame)
        self.scroll_area.setMinimumWidth(320)
        self.scroll_area.setMaximumWidth(450)

        sidebar_widget = QtWidgets.QWidget()
        self.sidebar_layout = QtWidgets.QVBoxLayout(sidebar_widget)
        self.scroll_area.setWidget(sidebar_widget)
        
        splitter.addWidget(self.scroll_area)

        # 3D Viewport Container
        self.plotter_widget = QtWidgets.QWidget()
        self.plotter_layout = QtWidgets.QVBoxLayout(self.plotter_widget)
        self.plotter_layout.setContentsMargins(0, 0, 0, 0)
        
        # Add the plotter itself (which is the QWidget) to the layout
        self.plotter = QtInteractor(self.plotter_widget)
        self.plotter_layout.addWidget(self.plotter)
        self.plotter.set_background("black")
        self.plotter.enable_eye_dome_lighting()
        
        splitter.addWidget(self.plotter_widget)
        
        # Ensure the 3D window gets the stretch
        splitter.setStretchFactor(1, 1)

        self.setup_sidebar(self.sidebar_layout)

    def setup_sidebar(self, layout):
        layout.addWidget(QtWidgets.QLabel("<b>FITS File</b>"))
        btn_open = QtWidgets.QPushButton("OPEN FITS...")
        btn_open.clicked.connect(self.open_file)
        layout.addWidget(btn_open)

        self.lbl_file = QtWidgets.QLabel("No file selected")
        self.lbl_file.setWordWrap(True)
        layout.addWidget(self.lbl_file)

        layout.addSpacing(10)

        # Preview Image
        self.lbl_preview = QtWidgets.QLabel()
        self.lbl_preview.setFixedSize(280, 200)
        self.lbl_preview.setAlignment(QtCore.Qt.AlignCenter)
        self.lbl_preview.setStyleSheet("border: 1px solid #555; background-color: #222;")
        self.lbl_preview.setText("No Preview")
        layout.addWidget(self.lbl_preview)

        layout.addSpacing(10)
        layout.addWidget(QtWidgets.QLabel("<b>Image Settings</b>"))
        
        self.spn_hdu = self.create_spinbox("HDU Index:", 0, 20, 0)
        self.spn_hdu[1].valueChanged.connect(self.update_preview)
        layout.addLayout(self.spn_hdu[0])

        self.spn_size = self.create_double_spinbox("Physical Size (mm):", 10.0, 500.0, 150.0)
        layout.addLayout(self.spn_size[0])

        self.spn_height = self.create_double_spinbox("Max Height (mm):", 1.0, 50.0, 10.0)
        layout.addLayout(self.spn_height[0])

        self.spn_base = self.create_double_spinbox("Base Thickness (mm):", 0.5, 20.0, 2.0)
        layout.addLayout(self.spn_base[0])

        self.spn_clip = self.create_double_spinbox("Clipping (%):", 0, 10.0, 1.0)
        self.spn_clip[1].valueChanged.connect(self.update_preview)
        layout.addLayout(self.spn_clip[0])
        
        self.spn_smooth = self.create_double_spinbox("Smoothing (sigma):", 0, 10.0, 1.5)
        self.spn_smooth[1].valueChanged.connect(self.update_preview)
        layout.addLayout(self.spn_smooth[0])
        
        self.spn_downsample = self.create_spinbox("Downsample Factor:", 1, 10, 2)
        self.spn_downsample[1].valueChanged.connect(self.update_preview)
        layout.addLayout(self.spn_downsample[0])

        self.cmb_scale = QtWidgets.QComboBox()
        self.cmb_scale.addItems(["log", "asinh", "linear"])
        self.cmb_scale.currentTextChanged.connect(self.update_preview)
        layout.addWidget(QtWidgets.QLabel("Scaling Mode:"))
        layout.addWidget(self.cmb_scale)
        
        self.chk_invert = QtWidgets.QCheckBox("Invert Heights")
        self.chk_invert.stateChanged.connect(self.update_preview)
        layout.addWidget(self.chk_invert)

        layout.addSpacing(10)
        
        # Border
        layout.addWidget(QtWidgets.QLabel("<b>Border</b>"))
        self.spn_border_w = self.create_double_spinbox("Width (mm):", 0, 20.0, 0.0)
        layout.addLayout(self.spn_border_w[0])
        self.spn_border_h = self.create_double_spinbox("Height (mm):", 0, 10.0, 0.0)
        layout.addLayout(self.spn_border_h[0])

        layout.addSpacing(10)
        layout.addWidget(QtWidgets.QLabel("<b>Custom Label</b>"))
        self.txt_label = QtWidgets.QLineEdit()
        self.txt_label.setPlaceholderText("Enter text for raised label...")
        layout.addWidget(self.txt_label)

        layout.addWidget(QtWidgets.QLabel("Label Style:"))
        self.cmb_label_mode = QtWidgets.QComboBox()
        self.cmb_label_mode.addItems(["Standard Text", "Braille"])
        layout.addWidget(self.cmb_label_mode)

        self.spn_label_depth = self.create_double_spinbox("Text Depth (mm):", 0.1, 5.0, 1.5)
        layout.addLayout(self.spn_label_depth[0])
        
        self.spn_label_size = self.create_double_spinbox("Font Size (mm):", 1.0, 50.0, 5.0)
        layout.addLayout(self.spn_label_size[0])

        self.spn_label_area = self.create_double_spinbox("Label Area (mm):", 0.0, 50.0, 10.0)
        layout.addLayout(self.spn_label_area[0])

        layout.addStretch()

        # Final Actions
        btn_recenter = QtWidgets.QPushButton("CENTER CAMERA")
        btn_recenter.clicked.connect(self.recenter_view)
        layout.addWidget(btn_recenter)

        self.btn_process = QtWidgets.QPushButton("PROCESS 3D MODEL")
        self.btn_process.setFixedHeight(50)
        self.btn_process.setStyleSheet("background-color: #3872a3; color: white; font-weight: bold;")
        self.btn_process.clicked.connect(self.start_processing)
        layout.addWidget(self.btn_process)

        self.progress_bar = QtWidgets.QProgressBar()
        self.progress_bar.setVisible(False)
        layout.addWidget(self.progress_bar)

        self.btn_save = QtWidgets.QPushButton("SAVE STL AS...")
        self.btn_save.clicked.connect(self.save_stl)
        self.btn_save.setEnabled(False)
        layout.addWidget(self.btn_save)

        self.status = QtWidgets.QLabel("Ready")
        layout.addWidget(self.status)

    def create_spinbox(self, label, min_v, max_v, default):
        lyt = QtWidgets.QHBoxLayout()
        lyt.addWidget(QtWidgets.QLabel(label))
        sb = QtWidgets.QSpinBox()
        sb.setRange(min_v, max_v); sb.setValue(default)
        lyt.addWidget(sb)
        return lyt, sb

    def create_double_spinbox(self, label, min_v, max_v, default):
        lyt = QtWidgets.QHBoxLayout()
        lyt.addWidget(QtWidgets.QLabel(label))
        sb = QtWidgets.QDoubleSpinBox()
        sb.setRange(min_v, max_v); sb.setSingleStep(0.1); sb.setValue(default)
        lyt.addWidget(sb)
        return lyt, sb

    def open_file(self):
        file_path, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Open FITS", "", "FITS (*.fits *.fit *.fz)")
        if file_path:
            self.fits_path = file_path
            self.lbl_file.setText(os.path.basename(file_path))
            self.raw_data = None  # Force reload
            self.update_preview()

    def update_preview(self):
        if not self.fits_path:
            return

        try:
            hdu_idx = self.spn_hdu[1].value()
            if self.raw_data is None or self.last_hdu != hdu_idx:
                self.raw_data = load_fits_data(self.fits_path, hdu_idx)
                self.last_hdu = hdu_idx

            # Preprocess image for preview
            norm_data = preprocess_image(
                self.raw_data,
                downsample_factor=self.spn_downsample[1].value(),
                clip_percentile=self.spn_clip[1].value(),
                log_scale=(self.cmb_scale.currentText() == 'log'),
                asinh_scale=(self.cmb_scale.currentText() == 'asinh'),
                invert=self.chk_invert.isChecked()
            )

            # Normalize to 0-255 for display
            img_8bit = (norm_data * 255).astype(np.uint8)
            h, w = img_8bit.shape
            
            # Create QImage from the numpy array
            qimg = QtGui.QImage(img_8bit.data, w, h, w, QtGui.QImage.Format_Grayscale8).copy()
            pixmap = QtGui.QPixmap.fromImage(qimg)
            
            # Scale to fit the preview label
            scaled_pixmap = pixmap.scaled(
                self.lbl_preview.size(), 
                QtCore.Qt.KeepAspectRatio, 
                QtCore.Qt.SmoothTransformation
            )
            self.lbl_preview.setPixmap(scaled_pixmap)
            
        except Exception as e:
            self.lbl_preview.setText(f"Preview Error: {str(e)}")

    def start_processing(self):
        if not self.fits_path: return
        params = {
            'fits_path': self.fits_path, 'hdu': self.spn_hdu[1].value(),
            'longest_side': self.spn_size[1].value(), 'max_height': self.spn_height[1].value(),
            'base_thickness': self.spn_base[1].value(), 'clip': self.spn_clip[1].value(),
            'smooth': self.spn_smooth[1].value(), 'downsample': self.spn_downsample[1].value(),
            'scale_mode': self.cmb_scale.currentText(), 'invert': self.chk_invert.isChecked(),
            'border_width': self.spn_border_w[1].value(), 'border_height': self.spn_border_h[1].value(),
            'label_text': self.txt_label.text().strip(),
            'label_mode': "braille" if self.cmb_label_mode.currentText() == "Braille" else "text",
            'label_depth': self.spn_label_depth[1].value(),
            'label_size': self.spn_label_size[1].value(),
            'label_area': self.spn_label_area[1].value()
        }
        self.btn_process.setEnabled(False); self.progress_bar.setVisible(True)
        self.progress_bar.setRange(0, 0)
        self.worker = Worker(params)
        self.worker.finished.connect(self.update_view)
        self.worker.error.connect(self.handle_error)
        self.worker.start()

    def handle_error(self, message):
        QtWidgets.QMessageBox.critical(self, "Error", message)
        self.btn_process.setEnabled(True); self.progress_bar.setVisible(False)

    def recenter_view(self, *args, mesh=None):
        """Forces the camera to look directly at the center from above."""
        # If called by button click, args might contain a boolean, so we only use explicit mesh kwarg
        # or we try to get it from the plotter's actors
        if mesh is None:
            try:
                actors = self.plotter.renderer.actors
                for actor in actors.values():
                    if hasattr(actor, 'mapper') and actor.mapper is not None:
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

        # Use a small delay (100ms) to ensure the camera fits correctly after any layout changes
        QtCore.QTimer.singleShot(100, do_recenter)

    def update_view(self, result):
        self.stl_path = result['path']
        self.btn_process.setEnabled(True); self.btn_save.setEnabled(True)
        self.progress_bar.setVisible(False)
        self.status.setText(f"Done! Faces: {result['faces']:,}")
        
        self.plotter.clear()
        mesh = pv.read(str(self.stl_path))
        
        # Center the mesh visually so camera resets work perfectly in all views
        # This only affects the preview, not the exported STL file
        mesh.translate(-np.array(mesh.center), inplace=True)
        
        # Add mesh to scene - Smooth shading and light gray
        self.plotter.add_mesh(mesh, color="lightgray", smooth_shading=True, name="model")
        
        # Reset camera using the precise bounds of the mesh
        self.recenter_view(mesh)

    def save_stl(self):
        if not self.stl_path: return
        save_path, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Save STL", "", "STL (*.stl)")
        if save_path:
            import shutil
            shutil.copy(self.stl_path, save_path)
            QtWidgets.QMessageBox.information(self, "Success", f"Saved to {save_path}")

if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    window = AstroTouchWindow()
    window.show()
    sys.exit(app.exec())
