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
SAMPLE_NUMBER = 17
FILE_NAME = "smarna.txt"

#SAMPLE_NUMBER = 5
#FILE_NAME = "test1.txt"
###


def sample(list_input, subsample_size):
    subsample = rnd.choices(list_input, k=subsample_size)
    return [i.split(" ") for i in subsample]


def generate_all_sxs(K):
    """Returns a set of all simplices in simplicial complex K."""
    sxi = set()
    for sx in K:
        if type(sx[0]) != tuple:
            sx = [tuple(i) for i in sx]
        n = len(sx)
        for i in range(n):
            sxi = sxi.union(combinations(sx, i + 1))
    return sxi


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

f_points = {}
for i,point in enumerate(tri.points):
    f_points[i] = (point[0], point[1])

f_h = {}
for i in f_points.items():
    f_h[i[0]] = f[i[1]]

# Iskanje kriticnih celice
cx_points = []
for i in triangles:
    cx_points.append((i[0], i[1], i[2]))
complex_K = utils.closure(cx_points)
critical_cells = utils.extend(complex_K, f_h)

# IZRIS GRAFA
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_trisurf(triang, Z=h, cmap=plt.cm.Spectral)

# iz tock v koordiante in visino
points = {}
for i in f_points.items():
    points[i[0]] = (i[1][0], i[1][1], f_h[i[0]])

for i in critical_cells[1]:
    c_x = []
    c_y = []
    c_h = []
    for j in list(i):
        c_x.append(float(points[j][0]))
        c_y.append(float(points[j][1]))
        c_h.append(float(points[j][2]))
    ax.scatter(c_x, c_y, c_h)

plt.show()

# RAZŠIRITEV FUNKCIJE VIŠINE NA SIMPLICIALNI KOMPLEKS

#  V funkcijah, ki jih bomo uporabljali, bodo morale biti točke podane kot numpy array
#  parov (tj. tuple-ov).
reshaped_points_x_y = np.asarray(reshaped_points_x_y, dtype=tuple)

#  V spremenljivko T shranimo triangulacijo, ki je enaka kot tri, le da so točke pari in
#  ne seznami dolžine 2.
T = np.ndarray.tolist(reshaped_points_x_y[tri.simplices])

#  Klic funkcije extract, ki vrne par V, C po opravljenih krajšanjih kritičnih simpleksih.
V1, C1 = extract(generate_all_sxs(T), f, 8000)
#print('C', len(C), 'C1', len(C1))

