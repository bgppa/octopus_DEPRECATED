#Specific solver for the 1-dim heat equation - add comments more
from math import pi, sin, exp

def alpha(j):
	j += 1
	return (2*pi*j)**2

def phi(j, x):
	j += 1
	return sin(2*pi*j*x)

def solver(a, basis_expansion, y, obs_number):
	# a, y are lists, the others are doubles
	time_observed = 0.01
	h = time_observed / (obs_number - 1)
	tmp_sum = 0
	i = 0
	j = 0
	for i in range(0, obs_number):
		tmp_sum = 0
		for j in range(0, basis_expansion):
			tmp_sum += a[j] * exp(-alpha(j)*time_observed)*phi(j, h*i)
		y[i] = tmp_sum


