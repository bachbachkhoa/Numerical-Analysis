# Author: Duynt

import numpy as np

def isDDM(A):
	'''Check diagonally dominant matrix. Return
	   0: Not diagonally dominant matrix
	   1: Row diagonal dominance
	   2: Column diagonal dominance'''
	d, rowCheck, colCheck = len(A), 1, 2
	for i in range(d):
		rowSum, colSum = 0, 0
		for j in range(d):
			rowSum += abs(A[i][j])
			colSum += abs(A[j][i])
		if 2*abs(A[i][i]) < rowSum:
			rowCheck += 1
		if 2*abs(A[i][i]) < colSum:
			colCheck += 1
	if rowCheck == 1 and colCheck == 2:
		return 1
	elif rowCheck != 1 and colCheck == 2:
		return 2
	elif rowCheck == 1 and colCheck != 2:
		return 1
	else:
		return 0

def GaussSeidel(a, b, x0, err, Nmax):
	# ax = b -> x = Bx + d
	check = isDDM(a)
	if check == 0:
		return 'Coefficient matrix is not diagonally dominant matrix'
	elif check == 1:
		n = len(a)
		B = np.zeros([n, n])
		for i in range(n):
			B[i, :] = -a[i, :]/a[i, i]
		np.fill_diagonal(B, np.zeros(n))
		err *= (1 - np.linalg.norm(B, np.Inf))/np.linalg.norm(B, np.Inf)
		k = 0
		while True:
			x1 = np.zeros([n, 1])
			for i in range(n):
				x1[i] = (b[i] - np.sum(np.dot(a[i, :i], x1[:i])) \
						- np.sum(np.dot(a[i, i + 1:], x0[i + 1:])))/a[i, i]
			k += 1
			if np.linalg.norm(x1 - x0, np.Inf) <= err or k >= Nmax:
				break
			else:
				x0 = x1

	return x1, k

# a = np.array([[10, 2, -1, 2], [1, 5, 1, 0], [1, -2, -5, 1], [3, 0, 0, 9]])
# b = np.array([[-4], [-1], [2], [10]])
# x0 = np.array([[0], [0], [0], [0]])
# x, k = GaussSeidel(a, b, x0, 1e-5, 150)
# print(x)
# print(k)
