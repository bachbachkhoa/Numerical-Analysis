# Author: Duynt

import numpy as np

def quadraticFunction(data):
	n = len(data)
	X = data[:, 0]
	Y = data[:, 1]
	x = [n, sum([data[i, 0] for i in range(n)]), sum([data[i, 0]**2 for i in range(n)])]
	y = [sum([data[i, 0] for i in range(n)]), sum([data[i, 0]**2 for i in range(n)])\
		, sum([data[i, 0]**3 for i in range(n)])]
	z = [sum([data[i, 0]**2 for i in range(n)]), sum([data[i, 0]**3 for i in range(n)])\
		, sum([data[i, 0]**4 for i in range(n)])]
	t = [sum([data[i, 1] for i in range(n)]), sum([data[i, 1]*data[i, 0] for i in range(n)])\
		, sum([data[i, 1]*data[i, 0]**2 for i in range(n)])]

	P = np.flip(np.linalg.solve(np.array([x, y, z]), np.array(t)), 0)
	Pvalue = np.polyval(P, data[:, 0])
	err = np.sqrt(n**-1*sum([(data[i, 1] - Pvalue[i])**2 for i in range(n)]))
	
	return np.linalg.solve(np.array([x, y, z]), np.array(t)), err

def exponentialFunction(data):
	# y = a*e^bx
	n = len(data)
	data[:, 1] = np.log(data[:, 1])
	x = [n, sum(data[:, 0])]
	y = [sum(data[:, 0]), sum(data[:, 0]**2)]
	z = [sum(data[:, 1]), sum(data[:, 1]*data[:, 0])]
	P = np.linalg.solve([x, y], z)

	return np.exp(P[0]), P[1]

def polynomial(data):
	# y = a*x^b
	n = len(data)
	data[:, 1] = np.log(data[:, 1])
	data[:, 0] = np.log(data[:, 0])
	x = [n, sum(data[:, 0])]
	y = [sum(data[:, 0]), sum(data[:, 0]**2)]
	z = [sum(data[:, 1]), sum(data[:, 1]*data[:, 0])]
	P = np.linalg.solve([x, y], z)

	return np.exp(P[0]), P[1]
