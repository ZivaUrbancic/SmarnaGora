import random as rnd
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import Delaunay
from matplotlib.tri import Triangulation
from itertools import combinations, chain
# from collections import defaultdict
# from random import choice
# from math import sqrt
from cancel_critical_sx import extract
# from expand_function import extract_raw
# from readwrite import read_dgvf_from_file, write_dgvf_into_file

## Constants
SAMPLE_NUMBER = 2000
FILE_NAME = "smarna.txt"
###

def sample(list_input, subsample_size):
    subsample = rnd.choices(list_input, k=subsample_size)
    return [i.split(" ") for i in subsample]


f = {}
reshaped_points_x_y = []
data = open(FILE_NAME, "r").readlines()
points_x_y_h = sample(data, SAMPLE_NUMBER)
for i in points_x_y_h:
    reshaped_points_x_y.append((float(i[0]), float(i[1])))
    f[(float(i[0]), float(i[1]))] = int(i[2])
points_x_y = [i[:2] for i in points_x_y_h]
h = np.asarray([i[2] for i in points_x_y_h], dtype='float64')

# TRIANGULACIJA
tri = Delaunay(np.asarray(points_x_y, dtype='float64'))
triangles = tri.simplices
triang = Triangulation([i[0] for i in points_x_y_h], [i[1] for i in points_x_y_h], triangles=triangles)

# IZRIS GRAFA
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_trisurf(triang, Z=h, cmap=plt.cm.Spectral)

plt.show()

# RAZŠIRITEV FUNKCIJE VIŠINE NA SIMPLICIALNI KOMPLEKS


def generate_all_sxs(K):
    sxi = set()
    for sx in K:
        if type(sx[0]) != tuple:
            sx = [tuple(i) for i in sx]
        n = len(sx)
        for i in range(n):
            sxi = sxi.union(combinations(sx, i + 1))
    return sxi


reshaped_points_x_y = np. asarray(reshaped_points_x_y, dtype=tuple)
T = np.ndarray.tolist(reshaped_points_x_y[tri.simplices])

print(extract(generate_all_sxs(T), f, 4))
