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

# Serve the temp directory for STL preview
app.add_static_files('/temp', str(TEMP_DIR))

class AstroTouchGUI:
    def __init__(self):
        self.fits_path = None
        self.stl_path = None
        self.processing = False
        self.model_view = None
        self.scene = None
        self.stats_label = None
        self.loading_spinner = None
        
        # Parameters
        self.params = {
            'hdu': 0,
            'longest_side': 100.0,
            'max_height': 10.0,
            'base_thickness': 2.0,
            'clip': 1.0,
            'smooth': 1.5,
            'downsample': 2, # Default to 2 for safety
            'scale_mode': 'log',
            'invert': False,
            'border_width': 0.0,
            'border_height': 0.0
        }

    async def handle_upload(self, e):
        """Save uploaded FITS to a temporary file."""
        self.fits_path = TEMP_DIR / e.file.name
        content = await e.file.read()
        with open(self.fits_path, 'wb') as f:
            f.write(content)
        ui.notify(f'Uploaded {e.file.name}')
        
        # Try to auto-detect HDUs and check size
        try:
            with fits.open(self.fits_path) as hdul:
                hdus = [i for i, h in enumerate(hdul) if h.data is not None and h.data.ndim == 2]
                if hdus:
                    self.params['hdu'] = hdus[0]
                    shape = hdul[hdus[0]].data.shape
                    ui.notify(f'Image size: {shape[1]}x{shape[0]}')
                    if max(shape) > 1000 and self.params['downsample'] == 1:
                        ui.notify('Large image detected. Setting downsample to 2 for performance.', type='warning')
                        self.params['downsample'] = 2
        except:
            pass

    async def process(self):
        if not self.fits_path:
            ui.notify('Please upload a FITS file first', type='negative')
            return

        self.processing = True
        if self.loading_spinner:
            self.loading_spinner.set_visibility(True)
        
        try:
            # Run the heavy lifting in a thread
            result = await run.io_bound(self._run_processing)
            self.stl_path = result['path']
            
            # Update stats
            stats_text = f"Vertices: {result['vertices']:,} | Faces: {result['faces']:,}"
            if self.stats_label:
                self.stats_label.set_text(stats_text)
            
            ui.notify('Success! Model generated.')
            self.update_preview()
        except Exception as e:
            ui.notify(f'Error: {str(e)}', type='negative')
        finally:
            self.processing = False
            if self.loading_spinner:
                self.loading_spinner.set_visibility(False)

    def _run_processing(self):
        """The actual processing logic."""
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
        
        return {
            'path': output_path,
            'vertices': len(stl_mesh.vectors) * 3,
            'faces': len(stl_mesh.vectors)
        }

    def update_preview(self):
        if not self.stl_path or not self.model_view:
            return
            
        url = f'/temp/{self.stl_path.name}?t={os.path.getmtime(self.stl_path)}'
        
        self.model_view.delete()
        with self.scene:
            self.model_view = self.scene.group()
            with self.model_view:
                self.scene.stl(url).scale(0.1).move(z=-1)
                self.scene.spot_light(distance=100, intensity=0.8).move(y=-10, z=10)

    def download(self):
        if self.stl_path:
            # Direct static link is more robust for large files
            ui.download(f'/temp/{self.stl_path.name}')
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
            
            ui.number('Downsample Factor', min=1, max=10, step=1).bind_value(gui.params, 'downsample').classes('w-full')
            ui.label('Increasing this greatly speeds up large images.').classes('text-xs text-slate-500')
            
            ui.select({'log': 'Logarithmic', 'asinh': 'Arcsinh', 'linear': 'Linear'}, label='Scaling').bind_value(gui.params, 'scale_mode').classes('w-full')
            ui.checkbox('Invert Heights').bind_value(gui.params, 'invert')
            
            ui.separator().classes('my-4')
            
            ui.label('Border').classes('font-bold text-slate-600')
            ui.number('Border Width (mm)', step=1).bind_value(gui.params, 'border_width').classes('w-full')
            ui.number('Border Height (mm)', step=0.1).bind_value(gui.params, 'border_height').classes('w-full')

        # Main viewport
        with ui.column().classes('flex-grow h-full bg-slate-900 rounded-lg relative'):
            with ui.scene(width='100%', height='100%').classes('bg-slate-900') as scene:
                gui.scene = scene
                gui.model_view = scene.group()
                # Initial camera position
                scene.move_camera(x=0, y=-10, z=10, duration=0)
            
            # Info overlay
            with ui.column().classes('absolute top-4 left-4 text-white pointer-events-none'):
                ui.label('3D Preview').classes('text-lg opacity-50')
                gui.stats_label = ui.label('').classes('text-xs font-mono opacity-70')
            
            # Loading overlay
            with ui.column().classes('absolute inset-0 flex items-center justify-center bg-black/50') as gui.loading_spinner:
                ui.spinner(size='lg', color='primary')
                ui.label('Processing 3D Mesh...').classes('text-white mt-4')
            gui.loading_spinner.set_visibility(False)

if __name__ in {"__main__", "__mp_main__"}:
    ui.run(title='AstroTouch GUI', port=8080)
