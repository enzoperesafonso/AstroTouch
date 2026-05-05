#!/usr/bin/env python3
# AstroTouch: Convert astronomical FITS images to 3D printable STL files.
# Designed primarily for astronomy outreach for the blind and visually impaired.
#
# Original by Enzo Peres Afonso 2025
# Refactored and Enhanced

import argparse
import logging
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits
from stl import mesh
import pyvista as pv
from scipy.ndimage import gaussian_filter, zoom

# Try to import tqdm for progress bars, but don't fail if not present
try:
    from tqdm import tqdm
except ImportError:
    tqdm = None

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(levelname)s: %(message)s'
)
logger = logging.getLogger(__name__)

# Standard Grade 1 Braille mapping (a-z, 0-9, space, UEB punctuation)
BRAILLE_PATTERNS = {
    'a': [1,0,0,0,0,0], 'b': [1,1,0,0,0,0], 'c': [1,0,0,1,0,0], 'd': [1,0,0,1,1,0],
    'e': [1,0,0,0,1,0], 'f': [1,1,0,1,0,0], 'g': [1,1,0,1,1,0], 'h': [1,1,0,0,1,0],
    'i': [0,1,0,1,0,0], 'j': [0,1,0,1,1,0], 'k': [1,0,1,0,0,0], 'l': [1,1,1,0,0,0],
    'm': [1,0,1,1,0,0], 'n': [1,0,1,1,1,0], 'o': [1,0,1,0,1,0], 'p': [1,1,1,1,0,0],
    'q': [1,1,1,1,1,0], 'r': [1,1,1,0,1,0], 's': [0,1,1,1,0,0], 't': [0,1,1,1,1,0],
    'u': [1,0,1,0,0,1], 'v': [1,1,1,0,0,1], 'w': [0,1,0,1,1,1], 'x': [1,0,1,1,0,1],
    'y': [1,0,1,1,1,1], 'z': [1,0,1,0,1,1], ' ': [0,0,0,0,0,0],
    '1': [1,0,0,0,0,0], '2': [1,1,0,0,0,0], '3': [1,0,0,1,0,0], '4': [1,0,0,1,1,0],
    '5': [1,0,0,0,1,0], '6': [1,1,0,1,0,0], '7': [1,1,0,1,1,0], '8': [1,1,0,0,1,0],
    '9': [0,1,0,1,0,0], '0': [0,1,0,1,1,0], 
    '#': [0,0,1,1,1,1], # Number sign
    '^': [0,0,0,0,0,1], # Capital sign (Dot 6)
    '.': [0,1,0,0,1,1], ',': [0,1,0,0,0,0], '-': [0,0,1,0,0,1],
    '!': [0,1,1,0,1,0], '?': [0,1,1,0,0,1], ':': [0,1,0,0,1,0], ';': [0,1,1,0,0,0]
}

def generate_braille_mesh(text, dot_radius=0.75, dot_height=0.8, dot_spacing=2.5, cell_spacing=6.0):
    """Generates a 3D mesh of Braille dots for the given text (UEB compliant)."""
    meshes = []
    
    # 6-dot grid offsets (x, y)
    dot_offsets = [
        (0, 2*dot_spacing), (0, dot_spacing), (0, 0),       # Dots 1, 2, 3
        (dot_spacing, 2*dot_spacing), (dot_spacing, dot_spacing), (dot_spacing, 0) # Dots 4, 5, 6
    ]
    
    def add_cell(pattern, x_pos):
        for i, bit in enumerate(pattern):
            if bit:
                dx, dy = dot_offsets[i]
                dot = pv.Sphere(radius=dot_radius, center=(x_pos + dx, dy, 0))
                dot.scale([1, 1, dot_height/dot_radius], inplace=True)
                dot = dot.clip(normal='-z', origin=(0, 0, 0)).fill_holes(1000)
                meshes.append(dot)

    current_x = 0
    in_number_sequence = False
    
    i = 0
    while i < len(text):
        char = text[i]
        
        # 1. Number Indicator
        if char.isdigit():
            if not in_number_sequence:
                add_cell(BRAILLE_PATTERNS['#'], current_x)
                current_x += cell_spacing
                in_number_sequence = True
        else:
            in_number_sequence = False

        # 2. Capitalization Indicator
        if char.isupper():
            # Check for sequence of capitals (all-caps word)
            caps_count = 0
            j = i
            while j < len(text) and text[j].isupper():
                caps_count += 1
                j += 1
            
            if caps_count > 1:
                # Double capital sign for all-caps word
                add_cell(BRAILLE_PATTERNS['^'], current_x)
                current_x += cell_spacing
                add_cell(BRAILLE_PATTERNS['^'], current_x)
                current_x += cell_spacing
                # Add all capitalized characters
                for k in range(caps_count):
                    c = text[i+k].lower()
                    if c in BRAILLE_PATTERNS:
                        add_cell(BRAILLE_PATTERNS[c], current_x)
                        current_x += cell_spacing
                i += caps_count
                continue
            else:
                # Single capital sign
                add_cell(BRAILLE_PATTERNS['^'], current_x)
                current_x += cell_spacing
                char = char.lower()

        # 3. Main Character
        if char in BRAILLE_PATTERNS:
            add_cell(BRAILLE_PATTERNS[char], current_x)
            current_x += cell_spacing
        elif char.lower() in BRAILLE_PATTERNS:
            add_cell(BRAILLE_PATTERNS[char.lower()], current_x)
            current_x += cell_spacing
        
        i += 1
        
    if not meshes:
        return None
        
    return meshes[0].merge(meshes[1:]) if len(meshes) > 1 else meshes[0]

def load_fits_data(filepath, hdu_index=None):
    """
    Loads image data from a FITS file.
    If hdu_index is None, it attempts to find the first 2D image HDU.
    """
    logger.info(f"Loading FITS file: {filepath}")
    try:
        with fits.open(filepath) as hdul:
            if hdu_index is not None:
                if hdu_index >= len(hdul):
                    raise ValueError(f"HDU index {hdu_index} out of range. Max index: {len(hdul)-1}")
                hdu = hdul[hdu_index]
                if hdu.data is None or hdu.data.ndim != 2:
                    raise ValueError(f"HDU {hdu_index} does not contain 2D image data.")
            else:
                # Auto-detect HDU
                hdu = None
                for i, h in enumerate(hdul):
                    if h.data is not None and h.data.ndim == 2:
                        logger.info(f"Auto-detected 2D image data in HDU {i}")
                        hdu = h
                        break
                if hdu is None:
                    raise ValueError("Could not find any 2D image data in the FITS file.")
            
            # Return copy to avoid issues with closed memmap
            return hdu.data.astype(np.float32).copy()
    except Exception as e:
        logger.error(f"Could not read FITS file: {e}")
        raise

def preprocess_image(
    image_data,
    downsample_factor=1,
    clip_percentile=1.0,
    nan_value=None,
    log_scale=False,
    asinh_scale=False,
    invert=False
):
    """
    Applies various processing steps to the image data.
    Returns normalized data in range [0, 1].
    """
    if downsample_factor > 1:
        logger.info(f"Downsampling by factor {downsample_factor}...")
        image_data = zoom(image_data, 1 / downsample_factor, order=1)
        if image_data.size == 0:
            raise ValueError("Downsampling resulted in empty image.")

    ny, nx = image_data.shape
    logger.info(f"Processing image dimensions: {nx}x{ny}")

    finite_mask = np.isfinite(image_data)
    if not np.all(finite_mask):
        num_non_finite = np.sum(~finite_mask)
        logger.warning(f"Found {num_non_finite} NaN/inf values.")
        if nan_value is not None:
            image_data[~finite_mask] = nan_value
        else:
            # Placeholder, will be replaced with min after clipping
            image_data[~finite_mask] = 0

    # Clipping
    if clip_percentile and 0 < clip_percentile < 50:
        logger.info(f"Clipping data to {clip_percentile}% - {100-clip_percentile}% percentile range.")
        data_for_percentile = image_data[finite_mask] if not np.all(finite_mask) else image_data
        if data_for_percentile.size > 0:
            min_val = np.percentile(data_for_percentile, clip_percentile)
            max_val = np.percentile(data_for_percentile, 100 - clip_percentile)
            image_data = np.clip(image_data, min_val, max_val)
            # Update NaNs if they weren't explicitly set
            if not np.all(finite_mask) and nan_value is None:
                image_data[~finite_mask] = min_val
        else:
            logger.warning("Cannot calculate percentiles, no finite data available.")
    elif not np.all(finite_mask) and nan_value is None:
        # If no clipping but has NaNs, use global min
        min_val = np.nanmin(image_data)
        image_data[~finite_mask] = min_val

    # Scaling
    if log_scale:
        logger.info("Applying log1p scaling...")
        min_data = np.min(image_data)
        if min_data < 0:
            image_data -= min_data
        image_data = np.log1p(image_data)
    elif asinh_scale:
        logger.info("Applying arcsinh scaling...")
        min_data = np.min(image_data)
        # Shift to 0 and scale a bit (heuristic)
        image_data = np.arcsinh(image_data - min_data)

    if invert:
        logger.info("Inverting data height...")
        image_data = np.max(image_data) - image_data

    # Normalization to [0, 1] range
    min_val, max_val = np.min(image_data), np.max(image_data)
    data_range = max_val - min_val
    if data_range == 0:
        logger.warning("Data range is zero. Model will be flat.")
        normalized_data = np.zeros_like(image_data)
    else:
        normalized_data = (image_data - min_val) / data_range

    return normalized_data

def generate_mesh(
    normalized_data,
    longest_side_mm=None,
    max_height_mm=10.0,
    base_thickness_mm=2.0,
    smoothing_sigma=0,
    border_width_mm=0.0,
    border_height_mm=0.0,
    label_text="",
    label_depth_mm=1.5,
    label_font_size_mm=5.0,
    label_mode="text",
    label_area_mm=0.0
):
    """
    Generates a 3D mesh from the normalized image data.
    """
    ny, nx = normalized_data.shape
    
    # Calculate scaling factor (mm per pixel)
    if longest_side_mm and longest_side_mm > 0:
        scale_factor = longest_side_mm / max(nx, ny)
        logger.info(f"Scale factor: {scale_factor:.4f} mm/pixel")
    else:
        scale_factor = 1.0
        logger.info("Using 1 mm/pixel scaling (default)")

    # Apply height scaling and base thickness
    z_data = base_thickness_mm + (normalized_data * max_height_mm)

    # Smoothing
    if smoothing_sigma > 0:
        logger.info(f"Applying Gaussian smoothing (sigma={smoothing_sigma} pixels)...")
        z_data = gaussian_filter(z_data, sigma=smoothing_sigma)

    # Add border
    if border_width_mm > 0:
        border_pixels = int(round(border_width_mm / scale_factor))
        if border_pixels > 0:
            logger.info(f"Adding border: {border_width_mm}mm ({border_pixels} pixels)")
            ny_new = ny + 2 * border_pixels
            nx_new = nx + 2 * border_pixels
            # Border height is relative to the base thickness
            z_bordered = np.full((ny_new, nx_new), base_thickness_mm + border_height_mm, dtype=np.float32)
            z_bordered[border_pixels:border_pixels+ny, border_pixels:border_pixels+nx] = z_data
            z_data = z_bordered
            ny, nx = ny_new, nx_new
        else:
            logger.warning("Border width too small to represent at this scale, skipping.")

    # Add extra flat area for the label at the bottom
    if label_area_mm > 0:
        label_area_pixels = int(round(label_area_mm / scale_factor))
        if label_area_pixels > 0:
            logger.info(f"Adding dedicated label area: {label_area_mm}mm ({label_area_pixels} pixels)")
            # Prepend flat rows to the start of z_data (corresponds to negative Y)
            ny_new = ny + label_area_pixels
            z_label_area = np.full((ny_new, nx), base_thickness_mm + border_height_mm, dtype=np.float32)
            z_label_area[label_area_pixels:, :] = z_data
            z_data = z_label_area
            ny = ny_new

    # Create grid and center it around (0,0)
    x = (np.arange(nx, dtype=np.float32) - (nx - 1) / 2.0) * scale_factor
    y = (np.arange(ny, dtype=np.float32) - (ny - 1) / 2.0) * scale_factor
    xx, yy = np.meshgrid(x, y)

    num_vertices = nx * ny
    vertices = np.zeros((num_vertices * 2, 3), dtype=np.float32)
    
    # Top surface vertices
    vertices[:num_vertices, 0] = xx.flatten()
    vertices[:num_vertices, 1] = yy.flatten()
    vertices[:num_vertices, 2] = z_data.flatten()
    
    # Bottom surface vertices
    vertices[num_vertices:, 0] = xx.flatten()
    vertices[num_vertices:, 1] = yy.flatten()
    vertices[num_vertices:, 2] = 0.0

    logger.info("Generating mesh faces (vectorized)...")
    
    # Indices for grid cells
    i, j = np.meshgrid(np.arange(nx - 1), np.arange(ny - 1))
    v00 = (j * nx + i).flatten()
    v10 = (j * nx + (i + 1)).flatten()
    v01 = ((j + 1) * nx + i).flatten()
    v11 = ((j + 1) * nx + (i + 1)).flatten()

    # Top surface faces (pointing up)
    top_faces_1 = np.column_stack([v00, v10, v01])
    top_faces_2 = np.column_stack([v10, v11, v01])
    top_faces = np.vstack([top_faces_1, top_faces_2])

    # Bottom surface faces (pointing down)
    v00_b = v00 + num_vertices
    v10_b = v10 + num_vertices
    v01_b = v01 + num_vertices
    v11_b = v11 + num_vertices
    
    bottom_faces_1 = np.column_stack([v00_b, v01_b, v10_b])
    bottom_faces_2 = np.column_stack([v10_b, v01_b, v11_b])
    bottom_faces = np.vstack([bottom_faces_1, bottom_faces_2])

    # Side faces
    side_faces_list = []
    
    # Y-edges (top and bottom in image coordinates)
    for i_edge in [0, ny - 1]:
        idx = np.arange(nx - 1)
        v_t0 = i_edge * nx + idx
        v_t1 = i_edge * nx + idx + 1
        v_b0 = v_t0 + num_vertices
        v_b1 = v_t1 + num_vertices
        
        if i_edge == 0:
            side_faces_list.append(np.column_stack([v_t0, v_b0, v_b1]))
            side_faces_list.append(np.column_stack([v_t0, v_b1, v_t1]))
        else:
            side_faces_list.append(np.column_stack([v_t0, v_b1, v_b0]))
            side_faces_list.append(np.column_stack([v_t0, v_t1, v_b1]))

    # X-edges (left and right in image coordinates)
    for j_edge in [0, nx - 1]:
        idx = np.arange(ny - 1)
        v_t0 = idx * nx + j_edge
        v_t1 = (idx + 1) * nx + j_edge
        v_b0 = v_t0 + num_vertices
        v_b1 = v_t1 + num_vertices
        
        if j_edge == 0:
            side_faces_list.append(np.column_stack([v_t0, v_b1, v_b0]))
            side_faces_list.append(np.column_stack([v_t0, v_t1, v_b1]))
        else:
            side_faces_list.append(np.column_stack([v_t0, v_b0, v_b1]))
            side_faces_list.append(np.column_stack([v_t0, v_b1, v_t1]))

    all_faces = np.vstack([top_faces, bottom_faces] + side_faces_list)
    
    logger.info(f"Mesh generated: {len(vertices)} vertices, {len(all_faces)} faces")
    
    stl_mesh = mesh.Mesh(np.zeros(all_faces.shape[0], dtype=mesh.Mesh.dtype))
    stl_mesh.vectors = vertices[all_faces]
    
    # Add custom raised text label if requested
    if label_text:
        try:
            logger.info(f"Adding custom label: '{label_text}' (Mode: {label_mode}, Size: {label_font_size_mm}mm, Depth: {label_depth_mm}mm)")
            # 1. Generate label mesh
            if label_mode == "braille":
                # Braille dots use depth as the dot height
                text_mesh = generate_braille_mesh(label_text, dot_height=label_depth_mm)
            else:
                text_mesh = pv.Text3D(label_text, depth=label_depth_mm)
            
            if text_mesh is None:
                raise ValueError("Generated label mesh is empty.")
            
            # 2. Scale and position
            model_width = nx * scale_factor
            
            # Get text dimensions
            b = text_mesh.bounds
            text_w = b[1] - b[0]
            text_h = b[3] - b[2]
            
            # Initial scale
            if label_mode == "braille":
                scale = 1.0 # Braille is usually kept at physical scale
            else:
                scale = label_font_size_mm / text_h
            
            # Check if too wide (max 90% of model width)
            max_text_w = model_width * 0.9
            if (text_w * scale) > max_text_w:
                scale = max_text_w / text_w
                logger.info(f"Label too wide, scaling down to fit: {scale:.4f}")
            
            if scale != 1.0:
                text_mesh.scale([scale, scale, 1.0], inplace=True)
                text_w *= scale
                text_h *= scale
                
            # Center horizontally, place near bottom edge
            bottom_edge_y = -(ny - 1) / 2.0 * scale_factor
            
            b = text_mesh.bounds
            tx = - (b[1] + b[0]) / 2.0
            
            # If we have a dedicated label area, center the text within it
            if label_area_mm > 0:
                # The label area starts at bottom_edge_y and ends at bottom_edge_y + label_area_mm
                ty = bottom_edge_y + (label_area_mm / 2.0) - (text_h / 2.0) - b[2]
            else:
                ty = bottom_edge_y + (text_h * 0.2) - b[2]
                
            tz = base_thickness_mm + border_height_mm
            
            text_mesh.translate([tx, ty, tz], inplace=True)
            
            # 3. Convert pyvista mesh to numpy-stl
            tri_text = text_mesh.triangulate()
            text_faces = tri_text.faces.reshape(-1, 4)[:, 1:]
            text_vertices = tri_text.points
            
            text_stl_data = np.zeros(text_faces.shape[0], dtype=mesh.Mesh.dtype)
            text_stl = mesh.Mesh(text_stl_data)
            text_stl.vectors = text_vertices[text_faces]
            
            # 4. Merge meshes
            stl_mesh = mesh.Mesh(np.concatenate([stl_mesh.data, text_stl.data]))
            logger.info(f"Label added successfully ({len(text_faces)} faces)")
            
        except Exception as e:
            logger.error(f"Failed to add label: {e}")

    return stl_mesh

def main():
    parser = argparse.ArgumentParser(
        description="AstroTouch: Convert astronomical FITS images to 3D printable STL files.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Files
    parser.add_argument("input", type=Path, help="Input FITS file")
    parser.add_argument("output", type=Path, nargs='?', help="Output STL file (default: input_name.stl)")
    
    # Selection
    parser.add_argument("--hdu", type=int, help="HDU index (default: auto-detect first 2D image)")
    
    # Physical dimensions
    parser.add_argument("--longest-side", type=float, metavar="MM", help="Length of the longest side in mm")
    parser.add_argument("--max-height", type=float, default=10.0, metavar="MM", help="Maximum feature height above base")
    parser.add_argument("--base-thickness", type=float, default=2.0, metavar="MM", help="Thickness of the base plate")
    
    # Processing
    parser.add_argument("--clip", type=float, default=1.0, metavar="PERCENT", help="Percentile to clip outliers (0 to disable)")
    parser.add_argument("--smooth", type=float, default=0, metavar="SIGMA", help="Gaussian smoothing sigma in pixels (e.g., 1.5)")
    parser.add_argument("--downsample", type=int, default=1, metavar="FACTOR", help="Downsample factor (e.g., 2 halves dimensions)")
    parser.add_argument("--nan-value", type=float, help="Value to replace NaNs with (default: minimum)")
    parser.add_argument("--label", type=str, default="", help="Custom text label to add to the bottom edge")
    parser.add_argument("--label-mode", type=str, choices=["text", "braille"], default="text", help="Mode for the label (standard text or braille)")
    parser.add_argument("--label-depth", type=float, default=1.5, help="Depth (thickness) of the raised text in mm")
    parser.add_argument("--label-size", type=float, default=5.0, help="Font size (height) of the text in mm")
    parser.add_argument("--label-area", type=float, default=0.0, help="Width of dedicated flat area at the bottom for the label (mm)")
    
    # Scaling options
    scale_group = parser.add_mutually_exclusive_group()
    scale_group.add_argument("--log", action="store_true", help="Apply log(1+x) scaling to pixels")
    scale_group.add_argument("--asinh", action="store_true", help="Apply arcsinh scaling (often better for astronomy)")
    
    parser.add_argument("--invert", action="store_true", help="Invert height map (brightest pixels become lowest)")
    
    # Border
    parser.add_argument("--border-width", type=float, default=0.0, metavar="MM", help="Border width in mm")
    parser.add_argument("--border-height", type=float, default=0.0, metavar="MM", help="Border height above base in mm")
    
    # Misc
    parser.add_argument("--verbose", action="store_true", help="Enable verbose debug output")

    args = parser.parse_args()

    if args.verbose:
        logger.setLevel(logging.DEBUG)

    # Determine output path if not provided
    if not args.output:
        args.output = args.input.with_suffix('.stl')
        logger.info(f"No output path provided, using: {args.output}")

    try:
        # 1. Load data
        data = load_fits_data(args.input, args.hdu)
        
        # 2. Preprocess data (scaling, clipping, etc.)
        norm_data = preprocess_image(
            data,
            downsample_factor=args.downsample,
            clip_percentile=args.clip,
            nan_value=args.nan_value,
            log_scale=args.log,
            asinh_scale=args.asinh,
            invert=args.invert
        )
        
        # 3. Generate the 3D mesh
        stl_mesh = generate_mesh(
            norm_data,
            longest_side_mm=args.longest_side,
            max_height_mm=args.max_height,
            base_thickness_mm=args.base_thickness,
            smoothing_sigma=args.smooth,
            border_width_mm=args.border_width,
            border_height_mm=args.border_height,
            label_text=args.label,
            label_mode=args.label_mode,
            label_depth_mm=args.label_depth,
            label_font_size_mm=args.label_size,
            label_area_mm=args.label_area
        )
        
        # 4. Save the STL file
        logger.info(f"Saving STL to: {args.output}")
        stl_mesh.save(args.output)
        logger.info("Success! 3D model generated.")

    except KeyboardInterrupt:
        logger.info("\nInterrupted by user.")
        sys.exit(0)
    except Exception as e:
        logger.error(f"Error: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
