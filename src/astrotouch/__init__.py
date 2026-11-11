"""
AstroTouch: turn astronomical images into 3D-printable tactile models.

    from astrotouch import pipeline
    builder = pipeline.ModelBuilder.from_file("image.fits")
    settings = pipeline.apply_template(pipeline.Settings(), plaque="standard", recipe="radio")
    builder.mesh(settings, builder.height_map(settings)).save("model.stl")
"""
__version__ = "0.2.0"
