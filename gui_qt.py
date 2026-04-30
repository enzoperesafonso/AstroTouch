#!/usr/bin/env python3
import sys
import os
import threading
from pathlib import Path
import tempfile

import numpy as np
from astropy.io import fits
from PySide6 import QtWidgets, QtCore, QtGui
import pyvista as pv
from pyvistaqt import QtInteractor

# Import core logic
from fits_to_stl import load_fits_data, preprocess_image, generate_mesh

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
                border_height_mm=self.params['border_height']
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
        self.setWindowTitle("AstroTouch Desktop - 3D Astronomy")
        self.resize(1200, 800)

        self.fits_path = None
        self.stl_path = None

        # Main Layout
        central_widget = QtWidgets.QWidget()
        self.setCentralWidget(central_widget)
        main_layout = QtWidgets.QHBoxLayout(central_widget)

        # Sidebar
        sidebar = QtWidgets.QVBoxLayout()
        sidebar_widget = QtWidgets.QWidget()
        sidebar_widget.setFixedWidth(300)
        sidebar_widget.setLayout(sidebar)
        main_layout.addWidget(sidebar_widget)

        # 3D Viewport
        self.plotter = QtInteractor(self)
        main_layout.addWidget(self.plotter.interactor)
        self.plotter.set_background("black")
        self.plotter.add_axes()
        self.plotter.enable_shadows()

        # UI Elements
        self.setup_sidebar(sidebar)

    def setup_sidebar(self, layout):
        # File Section
        layout.addWidget(QtWidgets.QLabel("<b>Input FITS</b>"))
        btn_open = QtWidgets.QPushButton("Open FITS File")
        btn_open.clicked.connect(self.open_file)
        layout.addWidget(btn_open)
        self.lbl_file = QtWidgets.QLabel("No file selected")
        self.lbl_file.setWordWrap(True)
        layout.addWidget(self.lbl_file)

        layout.addSpacing(20)

        # Parameters
        layout.addWidget(QtWidgets.QLabel("<b>Dimensions</b>"))
        self.spn_hdu = self.create_spinbox("HDU Index:", 0, 10, 0)
        layout.addLayout(self.spn_hdu[0])
        
        self.spn_size = self.create_spinbox("Longest Side (mm):", 1, 500, 100)
        layout.addLayout(self.spn_size[0])
        
        self.spn_height = self.create_double_spinbox("Max Height (mm):", 0.1, 50.0, 10.0)
        layout.addLayout(self.spn_height[0])
        
        self.spn_base = self.create_double_spinbox("Base Thickness (mm):", 0.1, 10.0, 2.0)
        layout.addLayout(self.spn_base[0])

        layout.addSpacing(20)

        layout.addWidget(QtWidgets.QLabel("<b>Processing</b>"))
        self.spn_clip = self.create_double_spinbox("Clipping (%):", 0, 10.0, 1.0)
        layout.addLayout(self.spn_clip[0])
        
        self.spn_smooth = self.create_double_spinbox("Smoothing (sigma):", 0, 10.0, 1.5)
        layout.addLayout(self.spn_smooth[0])
        
        self.spn_downsample = self.create_spinbox("Downsample Factor:", 1, 10, 2)
        layout.addLayout(self.spn_downsample[0])

        self.cmb_scale = QtWidgets.QComboBox()
        self.cmb_scale.addItems(["log", "asinh", "linear"])
        layout.addWidget(QtWidgets.QLabel("Scaling Mode:"))
        layout.addWidget(self.cmb_scale)

        self.chk_invert = QtWidgets.QCheckBox("Invert Heights")
        layout.addWidget(self.chk_invert)

        layout.addSpacing(20)

        # Border
        layout.addWidget(QtWidgets.QLabel("<b>Border</b>"))
        self.spn_border_w = self.create_double_spinbox("Width (mm):", 0, 20.0, 0.0)
        layout.addLayout(self.spn_border_w[0])
        self.spn_border_h = self.create_double_spinbox("Height (mm):", 0, 10.0, 0.0)
        layout.addLayout(self.spn_border_h[0])

        layout.addStretch()

        # Action Buttons
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

        # Status Bar
        self.status = QtWidgets.QLabel("Ready")
        layout.addWidget(self.status)

    def create_spinbox(self, label, min_v, max_v, default):
        lyt = QtWidgets.QHBoxLayout()
        lyt.addWidget(QtWidgets.QLabel(label))
        sb = QtWidgets.QSpinBox()
        sb.setRange(min_v, max_v)
        sb.setValue(default)
        lyt.addWidget(sb)
        return lyt, sb

    def create_double_spinbox(self, label, min_v, max_v, default):
        lyt = QtWidgets.QHBoxLayout()
        lyt.addWidget(QtWidgets.QLabel(label))
        sb = QtWidgets.QDoubleSpinBox()
        sb.setRange(min_v, max_v)
        sb.setSingleStep(0.1)
        sb.setValue(default)
        lyt.addWidget(sb)
        return lyt, sb

    def open_file(self):
        file_path, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Open FITS File", "", "FITS Files (*.fits *.fit *.fz)")
        if file_path:
            self.fits_path = file_path
            self.lbl_file.setText(os.path.basename(file_path))
            self.status.setText(f"Loaded {os.path.basename(file_path)}")

    def start_processing(self):
        if not self.fits_path:
            QtWidgets.QMessageBox.warning(self, "Error", "Please select a FITS file first.")
            return

        params = {
            'fits_path': self.fits_path,
            'hdu': self.spn_hdu[1].value(),
            'longest_side': self.spn_size[1].value(),
            'max_height': self.spn_height[1].value(),
            'base_thickness': self.spn_base[1].value(),
            'clip': self.spn_clip[1].value(),
            'smooth': self.spn_smooth[1].value(),
            'downsample': self.spn_downsample[1].value(),
            'scale_mode': self.cmb_scale.currentText(),
            'invert': self.chk_invert.isChecked(),
            'border_width': self.spn_border_w[1].value(),
            'border_height': self.spn_border_h[1].value()
        }

        self.btn_process.setEnabled(False)
        self.progress_bar.setVisible(True)
        self.progress_bar.setRange(0, 0) # Indeterminate
        
        self.worker = Worker(params)
        self.worker.progress.connect(self.status.setText)
        self.worker.error.connect(self.handle_error)
        self.worker.finished.connect(self.update_view)
        self.worker.start()

    def handle_error(self, message):
        QtWidgets.QMessageBox.critical(self, "Error", message)
        self.btn_process.setEnabled(True)
        self.progress_bar.setVisible(False)
        self.status.setText("Error occurred")

    def update_view(self, result):
        self.stl_path = result['path']
        self.btn_process.setEnabled(True)
        self.btn_save.setEnabled(True)
        self.progress_bar.setVisible(False)
        self.status.setText(f"Done! Faces: {result['faces']:,}")

        # Load into PyVista
        self.plotter.clear()
        mesh = pv.read(str(self.stl_path))
        
        # Robust centering: manually shift all points to (0, 0, 0)
        # We center X and Y, but keep the bottom of the base at Z=0
        cx, cy, cz = mesh.center
        bounds = mesh.bounds # [xmin, xmax, ymin, ymax, zmin, zmax]
        
        # Translation vector: center in X/Y, but set Z-min to 0
        translation = [-cx, -cy, -bounds[4]]
        mesh.points += translation
        
        # Add mesh with smooth shading
        self.plotter.add_mesh(mesh, color="lightgray", show_edges=False, smooth_shading=True, name="model")
        
        # Enable Eye Dome Lighting for depth perception
        self.plotter.enable_eye_dome_lighting()
        
        # Visual aids: Floor grid and axes
        self.plotter.add_floor_grid(color="gray", line_width=1)
        self.plotter.add_axes()
        
        # Reset camera to look at the new center
        self.plotter.camera.focal_point = (0, 0, (bounds[5]-bounds[4])/2)
        self.plotter.reset_camera()

    def save_stl(self):
        if not self.stl_path:
            return
        
        save_path, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Save STL", "", "STL Files (*.stl)")
        if save_path:
            import shutil
            shutil.copy(self.stl_path, save_path)
            QtWidgets.QMessageBox.information(self, "Success", f"File saved to {save_path}")

if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    window = AstroTouchWindow()
    window.show()
    sys.exit(app.exec())
