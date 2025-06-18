AstroTouch
![alt text](https://img.shields.io/badge/License-MIT-yellow.svg)

<!-- Add other badges if you set up CI/CD, etc. -->
![alt text](m51.png)

Convert astronomical FITS images into 3D printable STL surface relief models, designed primarily for astronomy outreach initiatives for the blind and visually impaired (BVI). This script allows users to transform the brightness variations in a 2D FITS image into height variations on a tangible 3D model.

Motivation
Astronomical data is overwhelmingly visual. This project aims to bridge that gap by providing a tool to create tactile representations of celestial objects and data, making astronomy more accessible and engaging for individuals who are blind or visually impaired. By feeling the contours, peaks, and valleys corresponding to stars, nebulae, and galaxies, users can gain a different understanding of astronomical structures.

Features
Reads standard 2D FITS image data (.fits, .fit).
Supports compressed FITS files (.fits.fz, .fit.fz) automatically via astropy.
Allows selection of the specific HDU containing the image data.
Maps pixel brightness values to surface height (Z-axis).
Options for data processing to enhance tactile feel and printability:
Logarithmic Scaling (--log_scale): Enhances faint features, crucial for most astronomical images.
Clipping (--clip): Removes extreme outlier pixel values (e.g., saturated stars, cosmic rays) before scaling, preventing them from dominating the height range.
Smoothing (--smooth): Applies Gaussian smoothing to reduce noise and sharp pixel edges for a better tactile feel and significantly improved printability (reduces steep overhangs).
Inversion (--invert): Maps bright pixels to low points (pits) instead of high points (peaks).
Downsampling (--downsample): Reduces image resolution to simplify the model, speed up processing, reduce file size, and improve printability for very large images.
Adds a solid, flat base underneath the surface for print stability.
Adjustable maximum feature height (--max_height) and base thickness (--base_thickness) in millimeters.
Handles NaN/inf values in FITS data by replacing them with a sensible default (usually the post-clipping minimum).
Outputs standard STL files compatible with 3D printing slicers (like Cura, PrusaSlicer, etc.).
Example Workflow:
Here's a workflow example using an observation of the Helix Nebula in Hydrogen-alpha:

![alt text](lco_h_alpha_preview.jpg)

1. Input FITS Image:

Object: Helix Nebula (NGC 7293)
Filter: Hydrogen-alpha (Hα) - This traces ionized hydrogen gas, highlighting specific structures within the nebula.
Source: Las Cumbres Observatory (LCO) Global Telescope Network
Observation: Taken by Roly Warry on Wed Jun 12 2024, using the 0.4-meter telescope at Siding Spring Observatory, Australia.
File: (Let's assume a filename like coj0m421-sq37-20240612-0129-e91.fits.fz)
2.  Generate STL using recommended options:

Generated bash
python3 fits_to_stl.py coj0m421-sq37-20240612-0129-e91.fits.fz output_model.stl --hdu 1 --log_scale --clip 1 --max_height 50 --smooth 2.0 --downsample 2
content_copy
download
Use code with caution.
Bash
4.  3D Print the Model:

![alt text](both.png)

Installation
Prerequisites
Python 3.7+
pip (Python package installer)
Git (for cloning the repository)
Steps
Clone the repository:
Generated bash
git clone https://github.com/enzoperesafonso/AstroTouch.git
cd AstroTouch
content_copy
download
Use code with caution.
Bash
Create a virtual environment (Recommended):
Generated bash
python -m venv venv
# Activate the environment:
# On Windows:
# venv\Scripts\activate
# On macOS/Linux:
# source venv/bin/activate
content_copy
download
Use code with caution.
Bash
Install dependencies:
The required libraries are listed in requirements.txt.
Generated bash
pip install -r requirements.txt
content_copy
download
Use code with caution.
Bash
Usage
Run the script from your terminal (make sure your virtual environment is active):

Generated bash
python fits_to_stl.py <input_fits_file> <output_stl_file> [options]
content_copy
download
Use code with caution.
Bash
Required Arguments
<input_fits_file>: Path to the input FITS file (e.g., image.fits or image.fits.fz).
<output_stl_file>: Path where the output STL file will be saved (e.g., model.stl).
Optional Arguments
See python fits_to_stl.py --help for a full list and defaults. Key options include:

--hdu INDEX (Default: 0, often need 1)
--max_height MM (Default: 10.0)
--base_thickness MM (Default: 2.0)
--invert
--log_scale
--clip PERCENT (Default: 1.0)
--smooth SIGMA (Default: 0)
--downsample FACTOR (Default: 1)
Command Line Examples
Basic conversion (likely needing HDU 1):
Generated bash
python fits_to_stl.py my_image.fits.fz my_model.stl --hdu 1
content_copy
download
Use code with caution.
Bash
Recommended starting point for good tactile feel & printability:
Generated bash
python fits_to_stl.py nebula.fits nebula_tactile.stl --hdu 1 --log_scale --clip 1.0 --smooth 1.5 --max_height 12.0 --base_thickness 2.0
content_copy
download
Use code with caution.
Bash
Tips for Effective Models
(See the detailed tips in other sections of this README and within the script's help message: python fits_to_stl.py --help)

For Tactile Feel: Use --log_scale, --smooth, adjust --max_height and --clip. Simplify complex images with --downsample.
For 3D Printing: Use --smooth (crucial!), --base_thickness >= 2.0, consider --downsample. In your slicer: use a Brim, no supports on the tactile surface, scale appropriately, use PLA.
Troubleshooting
(See the detailed troubleshooting steps in other sections of this README and within the script's help message)

Common issues include needing --hdu 1, adjusting --smooth for printability, using a brim for bed adhesion, and using --downsample for large files.
License
This project is licensed under the MIT License - see the LICENSE file for details.

Acknowledgements
This tool relies heavily on Astropy, NumPy, SciPy, and NumPy-STL.
Inspired by efforts to make science accessible to everyone.
![alt text](ngc1300.png)
