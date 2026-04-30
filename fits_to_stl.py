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
        sys.exit(1)

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
    border_height_mm=0.0
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

    # Create grid
    x = np.arange(nx, dtype=np.float32) * scale_factor
    y = np.arange(ny, dtype=np.float32) * scale_factor
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
            border_height_mm=args.border_height
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
