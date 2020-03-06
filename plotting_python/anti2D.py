# Given a file where each list contains two columns,
# plot the data interpreting each line as a 2d point
# with x and y coordinates given by the columns
# You can pass a second parameter to be used as title

import numpy as np
import matplotlib
from matplotlib import pyplot as plt
import sys

"""Get the data from a txt file, where each line ends with a newline"""
f = open(sys.argv[1])
x = []
for i in f:
    x.append([float(a) for a in i[0:len(i)-1].split(" ")])
f.close()

# x should be now a list of 2-dimensional data
xx = [ i[1] for i in x]
yy = [ i[0] for i in x] 

try:
    plt.title(sys.argv[2])
except:
    pass

plt.plot(xx, yy, '*')
plt.show()
