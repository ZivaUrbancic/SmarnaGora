# Mogoce nama prav pride tekom dela: https://matplotlib.org/mpl_toolkits/mplot3d/tutorial.html#mpl_toolkits.mplot3d.Axes3D.quiver

import random as rnd
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

def sample(list_input, subsample_size):
    subsample = rnd.choices(list_input, k = subsample_size)
    return [i.split(" ") for i in subsample]

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

file_lines = open("smarna.txt","r").readlines()
file_line_subsample = sample(file_lines, 1000)

x =[float(i[0]) for i in file_line_subsample]
y =[float(i[1]) for i in file_line_subsample]
z =[float(i[2]) for i in file_line_subsample]
ax.scatter(x, y, z, c=z)

ax.set_zlabel('Visina')

plt.show()
