import pyvista as pv
import numpy as np
sphere = pv.Sphere()
center1 = sphere.center
sphere.translate(-np.array(center1), inplace=True)
center2 = sphere.center
print("Center 1:", center1)
print("Center 2:", center2)
