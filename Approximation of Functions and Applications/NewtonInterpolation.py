# -*- coding: utf-8 -*-
from __future__ import unicode_literals
import numpy as np
from sympy import *

x = symbols('x')

def dividedDifferenceTable(data):

	'''Input: 2-D array consists of 2 columns, the value column of x,
	          and the value column of y corresponding to x
	   Output: divided diffrence table'''

	n = len(data) - 1
	X = data[:, 0]

	f = zeros(n + 1)
	f[:, 0] = data[:, 1]
	for j in range(0, n):
		for i in range(0, n - j):
			f[i, j + 1] = (f[i + 1, j] - f[i, j])/(X[i + j + 1] - X[i])

	return f

def findX0(x, data):
	# find start index
	X = data[:, 0]
	return np.where(abs(x - X) == min(abs(x - X)))[0][0]

def NewtonInterpolation(f, c, data):

	'''Input: divided diffrence table
	   Output: polynomial interpolation P, P(c)
	'''

	n = len(data) - 1
	X = data[:, 0]
	startIndex = findX0(c, data)
	if startIndex <= n - startIndex:
		# Newton forward
		P = f[startIndex, 0]
		for i in range(1, n - startIndex + 1):
			t = 1
			for j in range(startIndex, startIndex + i):
				t *= (x - X[j])
			P += t*f[startIndex, i]
	else:
		# Newton backward
		P = f[startIndex, 0]
		for i in range(1, startIndex):
			t = 1
			for j in range(startIndex, startIndex - i, -1):
				t *= (x - X[j])
			P += t*f[start_index - i, i]

	P = simplify(P)
	P_c = P.subs(x, c)
	return P, P_c, startIndex



def main():
	# Input processing
	inFile = input("\nEnter input data file name: ")
	data = np.loadtxt(inFile)
	X = data[:, 0]
	n = len(X) - 1
	while True:
		print("\nEnter the value to be calculated in [", X[0], ", ", X[n], "]:", end = " ")
		c = float(input())
		if c >= X[0] and c <= X[n]:
			break
	
	f = dividedDifferenceTable(data)
	P, P_c, startIndex = NewtonInterpolation(f, c, data)
	print('\nP = ', P)
	print('P(', c, ') = ', P_c)

	# Do you want to enter more data?
	answer = input('Do you want to enter more data? (y/n)')
	while answer == 'y' or answer == 'Y':
		dataX = float(input('Value of x: '))
		dataY = float(input('Value of y: '))
		data = np.vstack([data, [dataX, dataY]])
		f = dividedDifferenceTable(data)
		P += simplify(np.prod(x - data[startIndex:, 0])*f[startIndex, len(data) - startIndex - 1])
		P_c = P.subs(x, c)
		print('\nP = ', P)
		print('P(', c, ') = ', P_c)


if __name__ == '__main__':
	main()


# def NewtonInterpolation1(data, c):

# 	'''Input: 2-D array consists of 2 columns, the value column of x,
# 	          and the value column of y corresponding to x, the value 
# 	          of x are evenly spaced; c in [x0, xn]
# 	   Output: P(c)
# 	'''

# 	n = len(data) - 1
# 	X = data[:, 0]
# 	Y = data[:, 1]

# 	# Create diffrence table
# 	f = zeros(n + 1)
# 	f[:, 0] = Y
# 	for j in range(n):
# 		for i in range(n - j):
# 			f[i, j + 1] = f[i + 1, j] - f[i, j]

# 	# Calculate h
# 	h = abs(X[1] - X[0])

# 	if abs(c - X[0]) <= abs(c - X[n]):
# 		# Newton forward interpolation
# 		t = (c - X[0])/h
# 		P = f[0, 0]
# 		for i in range(1, n + 1):
# 			s = 1
# 			for j in range(i):
# 				s *= (t - j)
# 			P += f[0, i]*s/factorial(i)
# 	else:
# 		# Newton backward interpolation
# 		t = (c - X[n])/h
# 		P = f[n, 0]
# 		for i in range(1, n + 1):
# 			s = 1
# 			for j in range(i):
# 				s *= (t + j)
# 			P += f[n - i, i]*s/factorial(i)

# 	return P
