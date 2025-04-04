# AstroTouch
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- Add other badges if you set up CI/CD, etc. -->

Convert astronomical FITS images into 3D printable STL surface relief models, designed primarily for astronomy outreach initiatives for the blind and visually impaired (BVI). This script allows users to transform the brightness variations in a 2D FITS image into height variations on a tangible 3D model.

## Motivation

Astronomical data is overwhelmingly visual. This project aims to bridge that gap by providing a tool to create tactile representations of celestial objects and data, making astronomy more accessible and engaging for individuals who are blind or visually impaired. By feeling the contours, peaks, and valleys corresponding to stars, nebulae, and galaxies, users can gain a different understanding of astronomical structures.

## Features

*   Reads standard 2D FITS image data (`.fits`, `.fit`).
*   Supports compressed FITS files (`.fits.fz`, `.fit.fz`) automatically via `astropy`.
*   Allows selection of the specific HDU containing the image data.
*   Maps pixel brightness values to surface height (Z-axis).
*   Options for data processing to enhance tactile feel and printability:
    *   **Logarithmic Scaling (`--log_scale`):** Enhances faint features, crucial for most astronomical images.
    *   **Clipping (`--clip`):** Removes extreme outlier pixel values (e.g., saturated stars, cosmic rays) before scaling, preventing them from dominating the height range.
    *   **Smoothing (`--smooth`):** Applies Gaussian smoothing to reduce noise and sharp pixel edges for a better tactile feel and significantly improved printability (reduces steep overhangs).
    *   **Inversion (`--invert`):** Maps bright pixels to low points (pits) instead of high points (peaks).
    *   **Downsampling (`--downsample`):** Reduces image resolution to simplify the model, speed up processing, reduce file size, and improve printability for very large images.
*   Adds a solid, flat base underneath the surface for print stability.
*   Adjustable maximum feature height (`--max_height`) and base thickness (`--base_thickness`) in millimeters.
*   Handles `NaN`/`inf` values in FITS data by replacing them with a sensible default (usually the post-clipping minimum).
*   Outputs standard STL files compatible with 3D printing slicers (like Cura, PrusaSlicer, etc.).

## Examples

*(Consider adding screenshots/photos here!)*

**Example Workflow:**

1.  Input FITS Image (e.g., a galaxy)
    ```
    [Insert Screenshot of Original FITS Galaxy Image Here]
    ```
2.  Generate STL using recommended options:
    ```bash
    python fits_to_stl.py galaxy.fits galaxy_model.stl --hdu 1 --log_scale --smooth 1.5 --clip 1.0 --max_height 12.0
    ```
3.  Preview STL in Slicer:
    ```
    [Insert Screenshot of STL Preview in Slicer Here]
    ```
4.  3D Print the Model:
    ```
    [Insert Photo of 3D Printed Galaxy Model Here]
    ```

## Installation

### Prerequisites

*   Python 3.7+
*   `pip` (Python package installer)
*   Git (for cloning the repository)

### Steps

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/[your-username]/[your-repo-name].git
    cd [your-repo-name]
    ```
    *(Replace `[your-username]/[your-repo-name]` with your actual GitHub username and repository name)*

2.  **Create a virtual environment (Recommended):**
    ```bash
    python -m venv venv
    # Activate the environment:
    # On Windows:
    # venv\Scripts\activate
    # On macOS/Linux:
    # source venv/bin/activate
    ```

3.  **Install dependencies:**
    The required libraries are listed in `requirements.txt`.
    ```bash
    pip install -r requirements.txt
    ```

## Usage

Run the script from your terminal (make sure your virtual environment is active):

```bash
python fits_to_stl.py <input_fits_file> <output_stl_file> [options]
