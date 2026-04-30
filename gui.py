#!/usr/bin/env python3
import base64
import io
import os
import tempfile
from pathlib import Path

from nicegui import app, run, ui
from astropy.io import fits
import numpy as np

# Import our core logic
from fits_to_stl import load_fits_data, preprocess_image, generate_mesh

# Configuration
TEMP_DIR = Path(tempfile.gettempdir()) / 'astrotouch'
TEMP_DIR.mkdir(parents=True, exist_ok=True)

class AstroTouchGUI:
    def __init__(self):
        self.fits_path = None
        self.stl_path = None
        self.processing = False
        self.model_view = None
        
        # Parameters
        self.params = {
            'hdu': 0,
            'longest_side': 100.0,
            'max_height': 10.0,
            'base_thickness': 2.0,
            'clip': 1.0,
            'smooth': 1.5,
            'downsample': 1,
            'scale_mode': 'log',
            'invert': False,
            'border_width': 0.0,
            'border_height': 0.0
        }

    async def handle_upload(self, e):
        """Save uploaded FITS to a temporary file."""
        self.fits_path = TEMP_DIR / e.name
        with open(self.fits_path, 'wb') as f:
            f.write(e.content.read())
        ui.notify(f'Uploaded {e.name}')
        
        # Try to auto-detect HDUs
        try:
            with fits.open(self.fits_path) as hdul:
                hdus = [i for i, h in enumerate(hdul) if h.data is not None and h.data.ndim == 2]
                if hdus:
                    self.params['hdu'] = hdus[0]
                    ui.notify(f'Found 2D data in HDU {hdus[0]}')
        except:
            pass

    async def process(self):
        if not self.fits_path:
            ui.notify('Please upload a FITS file first', type='negative')
            return

        self.processing = True
        ui.notify('Processing... this may take a few seconds.')
        
        try:
            # Run the heavy lifting in a thread to keep UI responsive
            self.stl_path = await run.cpu_bound(self._run_processing)
            ui.notify('Success! Model generated.')
            self.update_preview()
        except Exception as e:
            ui.notify(f'Error: {str(e)}', type='negative')
        finally:
            self.processing = False

    def _run_processing(self):
        """The actual processing logic called via run.cpu_bound."""
        data = load_fits_data(self.fits_path, self.params['hdu'])
        
        norm_data = preprocess_image(
            data,
            downsample_factor=self.params['downsample'],
            clip_percentile=self.params['clip'],
            log_scale=(self.params['scale_mode'] == 'log'),
            asinh_scale=(self.params['scale_mode'] == 'asinh'),
            invert=self.params['invert']
        )
        
        stl_mesh = generate_mesh(
            norm_data,
            longest_side_mm=self.params['longest_side'],
            max_height_mm=self.params['max_height'],
            base_thickness_mm=self.params['base_thickness'],
            smoothing_sigma=self.params['smooth'],
            border_width_mm=self.params['border_width'],
            border_height_mm=self.params['border_height']
        )
        
        output_path = TEMP_DIR / f"{self.fits_path.stem}.stl"
        stl_mesh.save(str(output_path))
        return output_path

    def update_preview(self):
        if not self.stl_path or not self.model_view:
            return
            
        # Convert STL to data URL for Three.js
        with open(self.stl_path, 'rb') as f:
            content = f.read()
            base64_data = base64.b64encode(content).decode('utf-8')
            data_url = f'data:application/sla;base64,{base64_data}'
        
        self.model_view.clear()
        with self.model_view:
            # Scale and center logic
            # We scale by 0.1 for the preview viewport
            ui.scene.stl(data_url).scale(0.1).move(z=-1)
            ui.scene.spot_light(distance=100, intensity=0.8).move(y=-10, z=10)

    def download(self):
        if self.stl_path:
            ui.download(self.stl_path)
        else:
            ui.notify('No model generated yet', type='warning')

@ui.page('/')
def main_page():
    gui = AstroTouchGUI()
    
    ui.colors(primary='#3872a3', secondary='#5091cd', accent='#00d1b2')

    with ui.header().classes('items-center justify-between bg-slate-800 text-white'):
        ui.label('AstroTouch 3D').classes('text-2xl font-bold')
        with ui.row():
            ui.button('Process', on_click=gui.process).bind_enabled_from(gui, 'processing', backward=lambda x: not x)
            ui.button('Download STL', icon='download', on_click=gui.download).classes('bg-green-600')

    with ui.row().classes('w-full h-screen no-wrap p-4'):
        # Sidebar for parameters
        with ui.column().classes('w-1/4 p-4 bg-slate-100 rounded-lg shadow-inner overflow-y-auto h-full'):
            ui.label('Configuration').classes('text-xl font-bold mb-4')
            
            ui.upload(label='Upload FITS', on_upload=gui.handle_upload, auto_upload=True).classes('w-full mb-4')
            
            ui.number('HDU Index', format='%d').bind_value(gui.params, 'hdu').classes('w-full')
            ui.number('Longest Side (mm)', step=1).bind_value(gui.params, 'longest_side').classes('w-full')
            ui.number('Max Feature Height (mm)', step=0.1).bind_value(gui.params, 'max_height').classes('w-full')
            ui.number('Base Thickness (mm)', step=0.1).bind_value(gui.params, 'base_thickness').classes('w-full')
            
            ui.separator().classes('my-4')
            
            ui.label('Image Processing').classes('font-bold text-slate-600')
            ui.slider(min=0, max=5, step=0.1).bind_value(gui.params, 'clip')
            ui.label().bind_text_from(gui.params, 'clip', backward=lambda x: f'Clipping: {x}%')
            
            ui.slider(min=0, max=10, step=0.5).bind_value(gui.params, 'smooth')
            ui.label().bind_text_from(gui.params, 'smooth', backward=lambda x: f'Smoothing (sigma): {x}')
            
            ui.select({'log': 'Logarithmic', 'asinh': 'Arcsinh', 'linear': 'Linear'}, label='Scaling').bind_value(gui.params, 'scale_mode').classes('w-full')
            ui.checkbox('Invert Heights').bind_value(gui.params, 'invert')
            
            ui.separator().classes('my-4')
            
            ui.label('Border').classes('font-bold text-slate-600')
            ui.number('Border Width (mm)', step=1).bind_value(gui.params, 'border_width').classes('w-full')
            ui.number('Border Height (mm)', step=0.1).bind_value(gui.params, 'border_height').classes('w-full')

        # Main viewport
        with ui.column().classes('flex-grow h-full bg-slate-900 rounded-lg relative'):
            with ui.scene(width='100%', height='100%').classes('bg-slate-900') as scene:
                gui.model_view = scene.group()
                # Initial camera position
                scene.move_camera(x=0, y=-10, z=10, duration=0)
            
            ui.label('3D Preview').classes('absolute top-4 left-4 text-white text-lg opacity-50 pointer-events-none')

if __name__ in {"__main__", "__mp_main__"}:
    ui.run(title='AstroTouch GUI', port=8080)
