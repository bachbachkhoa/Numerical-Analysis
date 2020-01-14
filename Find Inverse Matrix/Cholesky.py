# Author: Duynt

import numpy as np

'''Cholesky method to find inverse matrix
   Input: symmetric matrix
   Output: inverse matrix
'''

def isSquare (a):
	return all (len (row) == len (a) for row in a)

def checkSymmetric(a, tol = 1e-8):
	return np.all(np.abs(a-a.T) < tol)

def Cholesky(A):
	
	# A is symmetric matrix
	# Decomposite A = L*L^T such that L is lower triangular matrix
	L = np.zeros_like(A)
	n = len(L)
	for i in range(n):
		for j in range(i+1):
			if i==j:
				val = A[i][i] - np.sum(np.square(L[i][:i]))
				L[i][i] = np.sqrt(val)
			else:
				L[i][j] = (A[i][j] - np.sum(L[i][:j]*L[j][:j]))/L[j][j]

	# calculate a = L^-1
	a = np.zeros_like(A)
	np.fill_diagonal(a, [1/L[i, i] for i in range(n)])

	for j in range(n - 1):
		a[j + 1, j] = -L[j + 1, j]/L[j, j]/L[j + 1, j + 1]
		for i in range(j + 2, n):
			s = [L[i, k]*a[k, j] for k in range(1, i)]
			a[i, j] = -(np.sum(s) + L[i, j]/L[j, j])/L[i, i]

	return np.dot(a.T, a)

# A = np.array(A).astype(complex)
