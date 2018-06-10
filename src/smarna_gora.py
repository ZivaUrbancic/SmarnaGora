import random as rnd
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
from numpy import infty
from scipy.spatial import Delaunay
from matplotlib.tri import Triangulation, TriAnalyzer, UniformTriRefiner
from itertools import combinations, chain
# from collections import defaultdict
# from random import choice
# from math import sqrt
from cancel_critical_sx import extract, morse_complex
# from expand_function import extract_raw
# from readwrite import read_dgvf_from_file, write_dgvf_into_file
import utils, math

## Constants
SAMPLE_NUMBER = 1000
FILE_NAME = "smarna.txt"
np.random.seed(10)
# rnd.seed(45)
#SAMPLE_NUMBER = 5
#FILE_NAME = "test1.txt"
###

def sample(list_input, subsample_size):
    subsample = rnd.choices(list_input, k=subsample_size)
    return [i.split(" ") for i in subsample]

def bias_sample(list_input, subsample_size, limit = 500):
    list_input = list(map(lambda x : [float(x[0]),float(x[1]), float(x[2])], [i.split(" ") for i in list_input]))
    list_input_height_sorted = sorted(list_input, key = lambda x : x[2], reverse = True)
    upper = list(filter(lambda x : x[2] > limit, list_input_height_sorted))
    lower = list(filter(lambda x : x[2] <= limit, list_input_height_sorted))

    upper_len = len(upper)
    lower_len = len(lower)

    upper_percent = [1/(upper_len) * .05] * upper_len
    lower_percent = [1/(lower_len) * .95] * lower_len
    
    choice = np.random.choice(range(len(list_input_height_sorted)), subsample_size, p = upper_percent + lower_percent)
    return list(map(lambda x : list_input_height_sorted[x], choice))


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
points_x_y_h = bias_sample(data, SAMPLE_NUMBER)
for i in points_x_y_h:
    reshaped_points_x_y.append((float(i[0]), float(i[1])))
    f[(float(i[0]), float(i[1]))] = int(i[2])
points_x_y = [i[:2] for i in points_x_y_h]
h = np.asarray([i[2] for i in points_x_y_h], dtype='float64')

# TRIANGULACIJA
tri = Delaunay(np.asarray(points_x_y, dtype='float64'))
triangles = tri.simplices
print(len(tri.simplices))
triang = Triangulation([i[0] for i in points_x_y_h], [i[1] for i in points_x_y_h], triangles=triangles)

f_points = {}
for i,point in enumerate(tri.points):
    f_points[i] = (point[0], point[1])

f_h = {}
for i in f_points.items():
    f_h[i[0]] = f[i[1]]

points = sorted(points_x_y_h, key = lambda x: x[2], reverse = True)
print(points[-1], points[0])

# Iskanje kriticnih celice
# cx_points = []
# for i in triangles:
#     cx_points.append((i[0], i[1], i[2]))
# complex_K = utils.closure(cx_points)
# critical_cells = utils.extend(complex_K, f_h)

# IZRIS GRAFA
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
p = ax.plot_trisurf(triang, Z=h, cmap=plt.cm.Spectral)
# fig.colorbar(p)

# # iz tock v koordiante in visino
# points = {}
# for i in f_points.items():
#     points[i[0]] = (i[1][0], i[1][1], f_h[i[0]])


# RAZŠIRITEV FUNKCIJE VIŠINE NA SIMPLICIALNI KOMPLEKS

#  V funkcijah, ki jih bomo uporabljali, bodo morale biti točke podane kot numpy array
#  parov (tj. tuple-ov).
reshaped_points_x_y = np.asarray(reshaped_points_x_y, dtype=tuple)

#  V spremenljivko T shranimo triangulacijo, ki je enaka kot tri, le da so točke pari in
#  ne seznami dolžine 2.
T = np.ndarray.tolist(reshaped_points_x_y[tri.simplices])

#  Klic funkcije extract, ki vrne par V, C po opravljenih krajšanjih kritičnih simpleksih.
V1, C1 = extract(generate_all_sxs(T), f, 50)

critical_cells_smarna = []

for i in C1:
    c_x = []
    c_y = []
    c_h = []
    critical_cell_smarna = []
    for j in i:
        if (float(f[j]) > 1): # 500
            c_x.append(float(j[0]))
            c_y.append(float(j[1]))
            c_h.append(float(f[j]))
            critical_cell_smarna.append((float(j[0]), float(j[1]), float(f[j])))
    critical_cell_smarna = sorted(critical_cell_smarna, key = lambda x: x[2], reverse = True)
    if critical_cell_smarna:
        critical_cells_smarna.append(critical_cell_smarna)
    if(len(c_x) == 3):
        c_x.append(c_x[0])
        c_y.append(c_y[0])
        c_h.append(c_h[0])
    if(len(i) == 3):
        ax.plot(c_x, c_y, c_h, "k")
    elif(len(i) == 2):
        ax.plot(c_x, c_y, c_h, "r")
        
saddles = []
peaks = []

for i in critical_cells_smarna:
    if len(i) == 1:
        print("Minimum", i, sep=":\t")
    elif len(i) == 2:
        print("Sedlo", i, sep=":\t")
        saddles.append(i)
    elif len(i) == 3:
        print("Vrh", i, sep=":\t")
        peaks.append(i)
        
        
critical_cells_smarna = []

for i in C1:
    c_x = []
    c_y = []
    critical_cell_smarna = []
    for j in i:
        if (float(f[j]) > 500):
            c_x.append(float(j[0]))
            c_y.append(float(j[1]))
            critical_cell_smarna.append((float(j[0]), float(j[1])))
    if critical_cell_smarna:
        critical_cells_smarna.append(critical_cell_smarna)

# print(critical_cells_smarna)

def dopolni_simpleks(simpleks):
    tocke = []
    for tocka in simpleks:
        tocke.append((float(tocka[0]), float(tocka[1]), float(f[tocka])))
    return tuple(tocke)

paths = utils.find_paths_between_critical_simplices(V1, C1)
# print(critical_cells_smarna)
sedlo = ((458678.67, 109533.65), (458691.17, 109633.65))
# sedlo = ((457966.17, 109558.65), (458141.17, 109646.15))
barva = ["r", "b"]
indeks_barva = 0
for ccx in [a for a in critical_cells_smarna if len(a) == 3]:
    for i,j in paths.items():
        if (i[0] == tuple(ccx) and i[1] == sedlo) or (i[1] == tuple(ccx) and i[0] == sedlo):
            print("Pot:", j)
            trenutna_pot = []
            for pot in j[::-1]:
                c_x = []
                c_y = []
                c_h = []
                for i,korak in enumerate(pot):
                    simpleksi = dopolni_simpleks(korak)
                    for sim in simpleksi[:-1]:
                        if sim not in trenutna_pot:
                            trenutna_pot.append(sim)
                            c_x.append(float(sim[0]))
                            c_y.append(float(sim[1]))
                            c_h.append(float(sim[2]))
                        # print(sim)
                # ax.plot(c_x, c_y, c_h, barva[indeks_barva])
            indeks_barva += 1
                    # print(i,dopolni_simpleks(korak))


# # all_paths = []

# # for j,i in paths.items():
# #         for k in i:
# #             current_path = []
# #             # print(j, k, sep=": ")
# #             # print("POT------")
# #             for l in k:
# #                 if l not in current_path:
# #                     current_path.append(l)
# #                     # print(l)
# #             all_paths.append(current_path)
# #             # print("Dolzina pot: ", len(current_path))

# # peaks = sorted(peaks, key = lambda x: x, reverse = True)[:2]
# # saddles = sorted(saddles, key = lambda x: x, reverse = True)
            
# # # print(peaks)
# # # print(saddles)

# # peaks_top = []
# # for i in peaks:
# #     peaks_top.append(max(i, key = lambda x: x[2]))

# # saddles_top = []

# # for i in saddles:
# #     saddles_top.append(max(i, key = lambda x: x[2]))

# # peaks_top = sorted(peaks_top, key = lambda x: x[2], reverse = True)[:2]
# # # print(peaks_top)
# # middle_top = math.sqrt((peaks_top[0][0] - peaks_top[1][0])**2 + (peaks_top[0][1] - peaks_top[1][1])**2 + (peaks_top[0][2] - peaks_top[1][2])**2 ) / 2.

# # best_saddle_distance = 10000000
# # best_saddle = []
# # for i in saddles_top:
# #     left_peak_distance = math.sqrt((i[0] - peaks_top[1][0])**2 + (i[1] - peaks_top[1][1])**2 + (i[2] - peaks_top[1][2])**2 ) / 2.
# #     right_peak_distance = math.sqrt((peaks_top[0][0] - i[0])**2 + (peaks_top[0][1] - i[1])**2 + (peaks_top[0][2] - i[2])**2 ) / 2.
# #     distance = math.fabs(left_peak_distance - right_peak_distance)
# #     if distance < best_saddle_distance:
# #         best_saddle = i
# #         best_saddle_distance = distance

# # print(best_saddle, distance)

# # ax.scatter([best_saddle[0]], [best_saddle[1]], [best_saddle[2]], color = "r", s = [100])

# plt.tricontourf(triang, h, cmap=plt.cm.Spectral)
# plt.colorbar()
# plt.tricontour(triang, h, colors='k')

# plt.triplot([i[0] for i in points_x_y], [i[1] for i in points_x_y], tri.simplices.copy(), color = 'k', linewidth=0.7)

plt.show()