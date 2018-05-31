import random as rnd
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import Delaunay
from matplotlib.tri import Triangulation

## Constants
SAMPLE_NUMBER = 2000
FILE_NAME = "smarna.txt"
###

def sample(list_input, subsample_size):
    subsample = rnd.choices(list_input, k = subsample_size)
    return [i.split(" ") for i in subsample]

data = open(FILE_NAME, "r").readlines()

points_x_y_h = sample(data, SAMPLE_NUMBER)
points_x_y = [i[:2] for i in points_x_y_h]
h = np.asarray([i[2] for i in points_x_y_h], dtype='float64')

# TRIANGULACIJA
tri = Delaunay(np.asarray(points_x_y, dtype='float64'))
triangles = tri.simplices
triang = Triangulation([i[0] for i in points_x_y_h], [i[1] for i in points_x_y_h], triangles = triangles)

# IZRIS GRAFA
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_trisurf(triang, Z = h, cmap = plt.cm.Spectral)

plt.show()