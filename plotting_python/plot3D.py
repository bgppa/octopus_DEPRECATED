# Given a file where each line contains three columns,
# plot the data assuming that each line corresponds to a 3d point
# Each columns is interpreted as the x, y, and z coordinate respectively
# You can pass a title as a textual second shell parameter

import numpy as np
import sys
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

"""Get the data from a txt file, where each line ends with a newline"""
f = open(sys.argv[1])
x = []
for i in f:
    x.append([float(a) for a in i[0:len(i)-1].split(" ")])
f.close()

# x should be now a list of 3-dimensional data
#xx = [i[0] for i in x]
#yy = [i[1] for i in x]
#zz = [i[2] for i in x]
# Inverting the order because of the new format having the
# frequencies as first values in each row
zz = [i[0] for i in x]
xx = [i[1] for i in x]
yy = [i[2] for i in x]

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(xx, yy, zz, c='r', marker='o')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('time')

try:
    plt.title(sys.argv[2])
except:
    pass

plt.show()
