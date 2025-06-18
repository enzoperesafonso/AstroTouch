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
    longest_side_mm=None,
    max_height_mm=10.0,
    base_thickness_mm=2.0,
    invert=False,
    log_scale=False,
    clip_percentile=1.0,
    smoothing_sigma=0,
    downsample_factor=1,
    nan_value=None,
    border_width_mm=0.0,
    border_height_mm=0.0
):
    """
    Converts a 2D astronomical image from a FITS file into a 3D printable
    STL surface relief map.

    Args:
        fits_filepath (str): Path to the input FITS file.
        stl_filepath (str): Path to save the output STL file.
        hdu_index (int): Index of the HDU containing the image data.
        longest_side_mm (float, optional): The desired length of the longest side
                                           of the STL's base in millimeters.
                                           If None, the model's X/Y size will
                                           equal its pixel dimensions (1 pixel = 1 unit).
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
        nan_value (float, optional): Value to replace NaN/inf values with.
                                     If None, they are replaced with the minimum
                                     value of the data *after* clipping.
        border_width_mm (float): Adds a border of this width in millimeters.
                                 0 to disable. ### MODIFIED ###
        border_height_mm (float): The height of the border, measured from the base.
                                  E.g., 0.0 makes a flat flange.
    """

    print(f"--- Loading FITS file: {fits_filepath} (HDU {hdu_index}) ---")
    try:
        with fits.open(fits_filepath) as hdul:
            if len(hdul) <= hdu_index:
                raise ValueError(f"HDU index {hdu_index} out of range. Max index: {len(hdul)-1}")
            hdu = hdul[hdu_index]
            if hdu.data is None or hdu.data.ndim != 2:
                 raise ValueError(f"HDU {hdu_index} does not contain 2D image data.")
            image_data = hdu.data.astype(np.float32)
    except FileNotFoundError:
        print(f"ERROR: FITS file not found at {fits_filepath}")
        return
    except Exception as e:
        print(f"ERROR: Could not read FITS file or HDU: {e}")
        return

    print(f"Original image dimensions: {image_data.shape}")

    if downsample_factor > 1:
        print(f"Downsampling by factor {downsample_factor}...")
        image_data = zoom(image_data, 1 / downsample_factor, order=1)
        print(f"New image dimensions: {image_data.shape}")
        if image_data.size == 0:
            print("ERROR: Downsampling resulted in empty image.")
            return

    ny, nx = image_data.shape
    if nx < 2 or ny < 2:
        print("ERROR: Image dimensions are too small (< 2 pixels).")
        return

    finite_mask = np.isfinite(image_data)
    has_non_finite = not np.all(finite_mask)
    if has_non_finite:
        num_non_finite = np.sum(~finite_mask)
        print(f"WARNING: Found {num_non_finite} NaN/inf values.")
        if nan_value is not None:
            print(f"Replacing non-finite values with specified value: {nan_value}")
            image_data[~finite_mask] = nan_value
        else:
            print("Non-finite values will be replaced with the data minimum after clipping.")
            image_data[~finite_mask] = 0

    if clip_percentile and 0 < clip_percentile < 50:
        print(f"Clipping data to {clip_percentile:.2f}% - {100-clip_percentile:.2f}% percentile range.")
        data_for_percentile = image_data[finite_mask] if has_non_finite else image_data
        if data_for_percentile.size > 0:
            min_val = np.percentile(data_for_percentile, clip_percentile)
            max_val = np.percentile(data_for_percentile, 100 - clip_percentile)
            image_data = np.clip(image_data, min_val, max_val)
            print(f"Data clipped to range: [{min_val:.4g}, {max_val:.4g}]")
            if has_non_finite and nan_value is None:
                 print(f"Replacing original NaNs/Infs with clipped minimum: {min_val:.4g}")
                 image_data[~finite_mask] = min_val
        else:
            print("WARNING: Cannot calculate percentiles, no finite data available.")

    elif has_non_finite and nan_value is None:
         data_for_min = image_data[finite_mask]
         replace_val = np.min(data_for_min) if data_for_min.size > 0 else 0
         print(f"Replacing original NaNs/Infs with global minimum: {replace_val:.4g}")
         image_data[~finite_mask] = replace_val

    if log_scale:
        print("Applying log1p scaling (log(1+x))...")
        min_data = np.min(image_data)
        if min_data < 0:
            print(f"  Shifting data by {-min_data:.4g} to make it non-negative before log scaling.")
            image_data -= min_data
        image_data = np.log1p(image_data)

    if invert:
        print("Inverting data height...")
        image_data = -image_data

    print("Normalizing data for Z-axis...")
    min_val, max_val = np.min(image_data), np.max(image_data)
    data_range = max_val - min_val

    if data_range == 0:
        print("WARNING: Data range is zero after processing. Surface will be flat.")
        z_data = np.full(image_data.shape, base_thickness_mm, dtype=np.float32)
    else:
        z_data = base_thickness_mm + ((image_data - min_val) / data_range) * max_height_mm
        print(f"Data scaled to Z range: [{base_thickness_mm:.2f} mm, {base_thickness_mm + max_height_mm:.2f} mm]")

    if smoothing_sigma and smoothing_sigma > 0:
        print(f"Applying Gaussian smoothing with sigma={smoothing_sigma}...")
        z_data = gaussian_filter(z_data, sigma=smoothing_sigma)
    # --- [END UNCHANGED SECTIONS 1] ---

    # --- ### MODIFIED LOGIC (So far it works...): CALCULATE SCALING FACTOR FIRST ### ---
    # We calculate the scale factor based on the *image content* dimensions (ny, nx).
    # This factor (mm per pixel) is then used to determine the border width in pixels.
    print("Calculating model scaling factor (mm/pixel)...")
    if longest_side_mm and longest_side_mm > 0:
        # Get dimensions of the image data *before* adding a border
        image_ny, image_nx = z_data.shape
        max_pixel_dim = max(image_nx, image_ny)
        scale_factor = longest_side_mm / max_pixel_dim
        print(f"Scaling content ({max_pixel_dim} pixels) to {longest_side_mm:.2f} mm. Scale factor: {scale_factor:.4f} mm/pixel.")
    else:
        scale_factor = 1.0
        print("Using 1-to-1 pixel-to-millimeter scaling.")

    # --- ### MODIFIED LOGIC: Add Border ### ---
    if border_width_mm > 0:
        # Convert border width from mm to pixels using the calculated scale factor
        border_width_pixels = int(round(border_width_mm / scale_factor))

        if border_width_pixels > 0:
            print(f"Adding a {border_width_mm:.2f} mm border ({border_width_pixels} pixels) with height {border_height_mm:.2f} mm.")

            border_z_value = base_thickness_mm + border_height_mm

            # Get original data dimensions
            orig_ny, orig_nx = z_data.shape

            # Create new larger array for bordered data
            ny_new = orig_ny + 2 * border_width_pixels
            nx_new = orig_nx + 2 * border_width_pixels
            z_data_bordered = np.full((ny_new, nx_new), border_z_value, dtype=np.float32)

            # Copy original z_data into the center
            start = border_width_pixels
            z_data_bordered[start:start + orig_ny, start:start + orig_nx] = z_data

            # Update z_data and dimensions for mesh generation
            z_data = z_data_bordered
            ny, nx = ny_new, nx_new # These are the final dimensions for the mesh grid
            print(f"Total dimensions with border: {z_data.shape}")
        else:
            print(f"WARNING: Border width {border_width_mm} mm is too small to represent at the current scale. Border skipped.")

    # --- ### MODIFIED LOGIC: Print final dimensions ### ---
    final_x_mm = nx * scale_factor
    final_y_mm = ny * scale_factor
    print(f"Final model base dimensions: {final_x_mm:.2f} mm x {final_y_mm:.2f} mm.")

    # --- Mesh Generation ---
    # This part now uses the final 'nx' and 'ny' (which may include the border)
    # and the pre-calculated 'scale_factor' to create the physical mesh.
    print("Generating mesh vertices...")
    x = np.arange(nx, dtype=np.float32) * scale_factor
    y = np.arange(ny, dtype=np.float32) * scale_factor
    xx, yy = np.meshgrid(x, y)

    # --- [UNCHANGED SECTIONS 2: Mesh faces, saving, etc.] ---
    num_vertices = nx * ny
    vertices = np.zeros((num_vertices * 2, 3), dtype=np.float32)
    vertices[:num_vertices, 0] = xx.flatten()
    vertices[:num_vertices, 1] = yy.flatten()
    vertices[:num_vertices, 2] = z_data.flatten()
    vertices[num_vertices:, 0] = xx.flatten()
    vertices[num_vertices:, 1] = yy.flatten()
    vertices[num_vertices:, 2] = 0.0

    print("Generating mesh faces...")
    num_faces_per_surface = 2 * (nx - 1) * (ny - 1)
    num_side_faces = 2 * (nx - 1) + 2 * (ny - 1)
    total_faces = (num_faces_per_surface * 2) + (num_side_faces * 2)
    faces = np.zeros((total_faces, 3), dtype=np.uint32)

    face_idx = 0
    base_offset = num_vertices
    for j in range(ny - 1):
        for i in range(nx - 1):
            v00, v10 = j * nx + i, j * nx + (i + 1)
            v01, v11 = (j + 1) * nx + i, (j + 1) * nx + (i + 1)
            faces[face_idx] = [v00, v10, v01]; face_idx += 1
            faces[face_idx] = [v10, v11, v01]; face_idx += 1
    for j in range(ny - 1):
        for i in range(nx - 1):
            v00, v10 = base_offset + j * nx + i, base_offset + j * nx + (i + 1)
            v01, v11 = base_offset + (j + 1) * nx + i, base_offset + (j + 1) * nx + (i + 1)
            faces[face_idx] = [v00, v01, v10]; face_idx += 1
            faces[face_idx] = [v10, v01, v11]; face_idx += 1
    for i in range(nx - 1):
        v_t0, v_t1 = i, i + 1
        v_b0, v_b1 = base_offset + i, base_offset + i + 1
        faces[face_idx] = [v_t0, v_b0, v_b1]; face_idx += 1
        faces[face_idx] = [v_t0, v_b1, v_t1]; face_idx += 1
        v_t0, v_t1 = (ny - 1) * nx + i, (ny - 1) * nx + i + 1
        v_b0, v_b1 = base_offset + (ny - 1) * nx + i, base_offset + (ny - 1) * nx + i + 1
        faces[face_idx] = [v_t0, v_b1, v_b0]; face_idx += 1
        faces[face_idx] = [v_t0, v_t1, v_b1]; face_idx += 1
    for j in range(ny - 1):
        v_t0, v_t1 = j * nx, (j + 1) * nx
        v_b0, v_b1 = base_offset + j * nx, base_offset + (j + 1) * nx
        faces[face_idx] = [v_t0, v_b1, v_b0]; face_idx += 1
        faces[face_idx] = [v_t0, v_t1, v_b1]; face_idx += 1
        v_t0, v_t1 = j * nx + (nx - 1), (j + 1) * nx + (nx - 1)
        v_b0, v_b1 = base_offset + j * nx + (nx - 1), base_offset + (j + 1) * nx + (nx - 1)
        faces[face_idx] = [v_t0, v_b0, v_b1]; face_idx += 1
        faces[face_idx] = [v_t0, v_b1, v_t1]; face_idx += 1

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
    # --- Command-line arguments ---
    parser.add_argument("fits_file", help="Path to the input FITS file.")
    parser.add_argument("stl_file", help="Path for the output STL file.")
    parser.add_argument("--hdu", type=int, default=0, help="Index of the HDU containing the 2D image data.")
    parser.add_argument(
        "--longest_side", type=float, default=None, metavar="MM",
        help="The desired length of the *image content* on its longest side (in mm). "
             "The final model will be larger if a border is added. If not set, "
             "defaults to 1-to-1 pixel-to-mm scaling."
    )
    parser.add_argument("--max_height", type=float, default=20.0, help="Maximum height of features above the base (mm).")
    parser.add_argument("--base_thickness", type=float, default=5.0, help="Thickness of the base plate below the features (mm).")
    parser.add_argument("--invert", action="store_true", help="Invert the height map (brightest pixels become lowest points).")
    parser.add_argument("--log_scale", action="store_true", help="Apply log1p (log(1+x)) scaling to pixel values before normalization.")
    parser.add_argument("--clip", type=float, default=1.0, metavar="PERCENT", help="Clip pixel values at the given low/high percentile (0-49.9). 0 disables clipping.")
    parser.add_argument("--smooth", type=float, default=0, metavar="SIGMA", help="Apply Gaussian smoothing with the given sigma (in pixels). 0 disables smoothing.")
    parser.add_argument("--downsample", type=int, default=1, metavar="FACTOR", help="Downsample image by this integer factor before processing. 1 means no downsampling.")
    parser.add_argument("--nan_value", type=float, default=None, metavar="VALUE", help="Value to replace NaN/inf pixels with. If not set, it defaults to the minimum value after clipping.")

    # ### MODIFIED: Command-line arguments for the border in mm ###
    parser.add_argument("--border_width_mm", type=float, default=0.0, metavar="MM", help="Add a border around the model with this width in millimeters. 0 disables the border.")
    parser.add_argument("--border_height", type=float, default=0.0, metavar="MM", help="Set the height of the border, measured from the base plate (in mm). For a flat flange, use 0.")

    args = parser.parse_args()

    # --- Argument Validation ---
    if args.longest_side is not None and args.longest_side <= 0:
        parser.error("--longest_side must be a positive number.")
    # ... (other validations are unchanged) ...
    if args.border_width_mm < 0:
        parser.error("--border_width_mm cannot be negative.")
    if args.border_height < 0:
        parser.error("--border_height cannot be negative.")

    # --- Call the function with updated arguments ---
    create_stl_from_fits(
        fits_filepath=args.fits_file,
        stl_filepath=args.stl_file,
        hdu_index=args.hdu,
        longest_side_mm=args.longest_side,
        max_height_mm=args.max_height,
        base_thickness_mm=args.base_thickness,
        invert=args.invert,
        log_scale=args.log_scale,
        clip_percentile=args.clip if args.clip > 0 else None,
        smoothing_sigma=args.smooth if args.smooth > 0 else None,
        downsample_factor=args.downsample,
        nan_value=args.nan_value,
        border_width_mm=args.border_width_mm,
        border_height_mm=args.border_height
    )
