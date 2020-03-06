import numpy as np
from math import sqrt
import matplotlib
from matplotlib import pyplot as plt
#from scipy.cluster.vq import kmeans
import sys
from numpy import histogram

"""Get the data from a txt file, where each line ends with a newline"""
f = open(sys.argv[1])
x = []
for i in f:
    x.append([float(a) for a in i[0:len(i)-1].split(" ")])
f.close()

#print (kmeans(x, int(sqrt(len(x)))))

q = histogram(x, bins=20)
for i in range(len(q[0])):
    print(q[1][i], q[0][i])
