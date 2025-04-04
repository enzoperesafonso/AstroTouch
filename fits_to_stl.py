# Enzo Peres Afonso 2025
import argparse
import numpy as np
from astropy.io import fits
from stl import mesh
from scipy.ndimage import gaussian_filter, zoom

def create_stl_from_fits(
    fits_filepath,
    stl_filepath,
    hdu_index=0,
    max_height_mm=10.0,
    base_thickness_mm=2.0,
    invert=False,
    log_scale=False,
    clip_percentile=1.0,
    smoothing_sigma=0,
    downsample_factor=1,
    nan_value=0.0 # Value to replace NaNs with AFTER potential scaling
):
    """
    Converts a 2D astronomical image from a FITS file into a 3D printable
    STL surface relief map.

    Args:
        fits_filepath (str): Path to the input FITS file.
        stl_filepath (str): Path to save the output STL file.
        hdu_index (int): Index of the HDU containing the image data.
        max_height_mm (float): Maximum height of the features above the base (in mm).
        base_thickness_mm (float): Thickness of the base plate (in mm).
        invert (bool): If True, maps higher pixel values to lower Z (pits).
        log_scale (bool): If True, applies log1p scaling to the data.
        clip_percentile (float): Percentile (0-100) to clip high/low values.
                                 E.g., 1.0 clips the top 1% and bottom 1%.
                                 Set to 0 or None to disable clipping.
        smoothing_sigma (float): Sigma for Gaussian smoothing (in pixels).
                                 Set to 0 or None to disable smoothing.
        downsample_factor (int): Factor to downsample the image by (e.g., 2 means
                                 halving dimensions). 1 means no downsampling.
        nan_value (float): Value used to replace NaN/inf values in the original data
                           *before* scaling. Defaults to the minimum value after clipping.
    """
    print(f"--- Loading FITS file: {fits_filepath} (HDU {hdu_index}) ---")
    try:
        with fits.open(fits_filepath) as hdul:
            if len(hdul) <= hdu_index:
                raise ValueError(f"HDU index {hdu_index} out of range. Max index: {len(hdul)-1}")
            hdu = hdul[hdu_index]
            if hdu.data is None or hdu.data.ndim != 2:
                 raise ValueError(f"HDU {hdu_index} does not contain 2D image data.")
            image_data = hdu.data.astype(np.float32) # Use float32 for processing
    except FileNotFoundError:
        print(f"ERROR: FITS file not found at {fits_filepath}")
        return
    except Exception as e:
        print(f"ERROR: Could not read FITS file or HDU: {e}")
        return

    print(f"Original image dimensions: {image_data.shape}")

    # --- Optional Downsampling ---
    if downsample_factor > 1:
        print(f"Downsampling by factor {downsample_factor}...")
        # Use order=1 for bilinear interpolation, better than nearest neighbor (order=0)
        image_data = zoom(image_data, 1 / downsample_factor, order=1)
        print(f"New image dimensions: {image_data.shape}")
        if image_data.size == 0:
            print("ERROR: Downsampling resulted in empty image.")
            return

    ny, nx = image_data.shape
    if nx < 2 or ny < 2:
        print("ERROR: Image dimensions too small after downsampling (< 2 pixels).")
        return

    # --- Handle NaNs and Infs ---
    # Replace NaNs/Infs *before* scaling/clipping
    finite_mask = np.isfinite(image_data)
    if not np.all(finite_mask):
        num_non_finite = np.sum(~finite_mask)
        print(f"WARNING: Found {num_non_finite} NaN/inf values.")
        # Decide replacement value: use specified nan_value, or estimate if needed.
        # A common strategy is to replace with the minimum *finite* value found *later*
        replace_val = nan_value # Placeholder, will update after potential clipping
        print(f"Replacing NaNs/Infs with placeholder value for now.")
        image_data[~finite_mask] = 0 # Temporary replacement
    else:
        print("Data contains no NaNs or Infs.")


    # --- Optional Clipping ---
    if clip_percentile and 0 < clip_percentile < 50:
        print(f"Clipping data to {clip_percentile:.2f}% - {100-clip_percentile:.2f}% percentile range.")
        # Calculate percentiles on finite data only if there were NaNs
        data_for_percentile = image_data[finite_mask] if not np.all(finite_mask) else image_data
        if data_for_percentile.size > 0:
            min_val = np.percentile(data_for_percentile, clip_percentile)
            max_val = np.percentile(data_for_percentile, 100 - clip_percentile)
            image_data = np.clip(image_data, min_val, max_val)
            print(f"Data clipped to range: [{min_val:.4g}, {max_val:.4g}]")
            # Now set the final replacement value for NaNs if they existed
            if not np.all(finite_mask):
                 replace_val = min_val # Use the clipped minimum
                 image_data[~finite_mask] = replace_val
                 print(f"Final replacement value for NaNs/Infs: {replace_val:.4g}")
        else:
            print("WARNING: Cannot calculate percentiles, no finite data available after initial NaN check?")
    elif not np.all(finite_mask):
        # No clipping, but NaNs existed. Replace them with the global finite minimum.
         data_for_minmax = image_data[finite_mask]
         if data_for_minmax.size > 0:
             replace_val = np.min(data_for_minmax)
         else:
             replace_val = 0 # Fallback if somehow all data was non-finite
         image_data[~finite_mask] = replace_val
         print(f"Final replacement value for NaNs/Infs: {replace_val:.4g}")


    # --- Optional Log Scaling ---
    if log_scale:
        print("Applying log1p scaling (log(1+x))...")
        # Ensure data is non-negative for log scaling
        min_data = np.min(image_data)
        if min_data < 0:
            print(f"  Shifting data by {-min_data:.4g} to make it non-negative before log scaling.")
            image_data -= min_data # Shift data so min is 0
        image_data = np.log1p(image_data)
        print("Log scaling applied.")

    # --- Optional Inversion ---
    if invert:
        print("Inverting data height...")
        image_data = -image_data # Higher original values are now lower

    # --- Normalization ---
    print("Normalizing data...")
    min_val = np.min(image_data)
    max_val = np.max(image_data)
    data_range = max_val - min_val

    if data_range == 0:
        print("WARNING: Data range is zero after processing. Surface will be flat.")
        # Assign a flat height, respecting the base thickness
        z_data = np.full(image_data.shape, base_thickness_mm, dtype=np.float32)
    else:
        # Scale data to [0, max_height_mm] and add base thickness
        z_data = base_thickness_mm + ((image_data - min_val) / data_range) * max_height_mm
        print(f"Data scaled to Z range: [{base_thickness_mm:.2f} mm, {base_thickness_mm + max_height_mm:.2f} mm]")


    # --- Optional Smoothing ---
    if smoothing_sigma and smoothing_sigma > 0:
        print(f"Applying Gaussian smoothing with sigma={smoothing_sigma}...")
        z_data = gaussian_filter(z_data, sigma=smoothing_sigma)
        print("Smoothing applied.")


    # --- Mesh Generation ---
    print("Generating mesh vertices...")
    x = np.arange(nx)
    y = np.arange(ny)
    xx, yy = np.meshgrid(x, y)

    # Vertices: Top surface first, then bottom surface
    # Top surface: (x, y, z_data(x,y))
    # Bottom surface: (x, y, 0) - a flat base
    num_vertices = nx * ny
    vertices = np.zeros((num_vertices * 2, 3), dtype=np.float32)

    # Top surface vertices (indices 0 to num_vertices-1)
    vertices[:num_vertices, 0] = xx.flatten()
    vertices[:num_vertices, 1] = yy.flatten()
    vertices[:num_vertices, 2] = z_data.flatten()

    # Bottom surface vertices (indices num_vertices to 2*num_vertices-1)
    vertices[num_vertices:, 0] = xx.flatten()
    vertices[num_vertices:, 1] = yy.flatten()
    vertices[num_vertices:, 2] = 0.0 # Base is at z=0


    print("Generating mesh faces...")
    num_faces_per_surface = 2 * (nx - 1) * (ny - 1)
    num_side_faces = 2 * (nx - 1) * 2 + 2 * (ny - 1) * 2 # Front/Back + Left/Right sides (2 triangles per quad)
    total_faces = num_faces_per_surface * 2 + num_side_faces
    faces = np.zeros((total_faces, 3), dtype=np.uint32) # Use uint32 for indices

    face_idx = 0

    # Top surface faces
    for j in range(ny - 1):
        for i in range(nx - 1):
            # Vertex indices for the current quad
            v00 = j * nx + i       # Top-left
            v10 = j * nx + (i + 1) # Top-right
            v01 = (j + 1) * nx + i # Bottom-left
            v11 = (j + 1) * nx + (i + 1) # Bottom-right

            # Triangle 1 (top-left, top-right, bottom-left)
            faces[face_idx] = [v00, v10, v01]
            face_idx += 1
            # Triangle 2 (top-right, bottom-right, bottom-left)
            faces[face_idx] = [v10, v11, v01]
            face_idx += 1

    # Bottom surface faces (indices offset by num_vertices)
    # Note: Winding order is reversed to face downwards (optional, matters for some viewers/slicers)
    base_offset = num_vertices
    for j in range(ny - 1):
        for i in range(nx - 1):
            v00 = base_offset + j * nx + i
            v10 = base_offset + j * nx + (i + 1)
            v01 = base_offset + (j + 1) * nx + i
            v11 = base_offset + (j + 1) * nx + (i + 1)

            # Triangle 1 (reversed)
            faces[face_idx] = [v00, v01, v10]
            face_idx += 1
            # Triangle 2 (reversed)
            faces[face_idx] = [v10, v01, v11]
            face_idx += 1

    # Side faces
    # Front edge (y=0)
    for i in range(nx - 1):
        v_top_0 = i
        v_top_1 = i + 1
        v_bot_0 = base_offset + i
        v_bot_1 = base_offset + i + 1
        faces[face_idx] = [v_top_0, v_bot_0, v_bot_1] # Triangle 1
        face_idx += 1
        faces[face_idx] = [v_top_0, v_bot_1, v_top_1] # Triangle 2
        face_idx += 1

    # Back edge (y=ny-1)
    y_back = ny - 1
    for i in range(nx - 1):
        v_top_0 = y_back * nx + i
        v_top_1 = y_back * nx + i + 1
        v_bot_0 = base_offset + y_back * nx + i
        v_bot_1 = base_offset + y_back * nx + i + 1
        faces[face_idx] = [v_top_0, v_bot_1, v_bot_0] # Triangle 1 (reversed winding for outside face)
        face_idx += 1
        faces[face_idx] = [v_top_0, v_top_1, v_bot_1] # Triangle 2 (reversed winding)
        face_idx += 1

    # Left edge (x=0)
    for j in range(ny - 1):
        v_top_0 = j * nx
        v_top_1 = (j + 1) * nx
        v_bot_0 = base_offset + j * nx
        v_bot_1 = base_offset + (j + 1) * nx
        faces[face_idx] = [v_top_0, v_bot_1, v_bot_0] # Triangle 1 (reversed winding)
        face_idx += 1
        faces[face_idx] = [v_top_0, v_top_1, v_bot_1] # Triangle 2 (reversed winding)
        face_idx += 1

    # Right edge (x=nx-1)
    x_right = nx - 1
    for j in range(ny - 1):
        v_top_0 = j * nx + x_right
        v_top_1 = (j + 1) * nx + x_right
        v_bot_0 = base_offset + j * nx + x_right
        v_bot_1 = base_offset + (j + 1) * nx + x_right
        faces[face_idx] = [v_top_0, v_bot_0, v_bot_1] # Triangle 1
        face_idx += 1
        faces[face_idx] = [v_top_0, v_bot_1, v_top_1] # Triangle 2
        face_idx += 1

    if face_idx != total_faces:
         print(f"WARNING: Face count mismatch! Expected {total_faces}, generated {face_idx}")

    # --- Create and Save STL ---
    print(f"Creating STL mesh object ({len(vertices)} vertices, {len(faces)} faces)...")
    surface_mesh = mesh.Mesh(np.zeros(faces.shape[0], dtype=mesh.Mesh.dtype))
    surface_mesh.vectors = vertices[faces]

    print(f"Saving STL file to: {stl_filepath}")
    try:
        surface_mesh.save(stl_filepath)
        print("--- STL file saved successfully! ---")
    except Exception as e:
        print(f"ERROR: Could not save STL file: {e}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert a 2D FITS image into a 3D printable STL surface relief map.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("fits_file", help="Path to the input FITS file.")
    parser.add_argument("stl_file", help="Path for the output STL file.")
    parser.add_argument("--hdu", type=int, default=0, help="Index of the HDU containing the 2D image data.")
    parser.add_argument("--max_height", type=float, default=10.0, help="Maximum height of features above the base (mm).")
    parser.add_argument("--base_thickness", type=float, default=2.0, help="Thickness of the base plate below the features (mm).")
    parser.add_argument("--invert", action="store_true", help="Invert the height map (brightest pixels become lowest points).")
    parser.add_argument("--log_scale", action="store_true", help="Apply log1p (log(1+x)) scaling to pixel values before normalization.")
    parser.add_argument("--clip", type=float, default=1.0, metavar="PERCENT", help="Clip pixel values at the given low/high percentile before scaling (0-50). 0 disables clipping.")
    parser.add_argument("--smooth", type=float, default=0, metavar="SIGMA", help="Apply Gaussian smoothing with the given sigma (in pixels). 0 disables smoothing.")
    parser.add_argument("--downsample", type=int, default=1, metavar="FACTOR", help="Downsample image by this integer factor before processing. 1 means no downsampling.")
    parser.add_argument("--nan_value", type=float, default=None, metavar="VALUE", help="Value to replace NaN/inf pixels with. Defaults to the minimum value after clipping (or global minimum if no clipping).")

    args = parser.parse_args()

    # Validate arguments
    if args.max_height <= 0:
        parser.error("--max_height must be positive.")
    if args.base_thickness < 0:
        parser.error("--base_thickness cannot be negative.")
    if args.clip < 0 or args.clip >= 50:
        parser.error("--clip percentile must be between 0 (exclusive) and 50 (inclusive). Use 0 to disable.")
    if args.smooth < 0:
        parser.error("--smooth sigma cannot be negative.")
    if args.downsample < 1:
        parser.error("--downsample factor must be 1 or greater.")

    nan_replace_strategy = args.nan_value if args.nan_value is not None else "auto" # Use 'auto' internally if not specified

    create_stl_from_fits(
        fits_filepath=args.fits_file,
        stl_filepath=args.stl_file,
        hdu_index=args.hdu,
        max_height_mm=args.max_height,
        base_thickness_mm=args.base_thickness,
        invert=args.invert,
        log_scale=args.log_scale,
        clip_percentile=args.clip if args.clip > 0 else None, # Pass None if 0
        smoothing_sigma=args.smooth if args.smooth > 0 else None, # Pass None if 0
        downsample_factor=args.downsample,
        nan_value = nan_replace_strategy # Let the function handle 'auto' logic
    )
