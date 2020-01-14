#Author: Duynt

import numpy as np

def GaussJordan(A):
	n = len(A)
	I = np.eye(n)
	x = np.arange(n)
	for i in x:
		pivot, u, v = A[i, i], np.copy(A[i, :]), np.copy(I[i, :])
		A[i, :] /= pivot
		I[i, :] /= pivot
		y = np.delete(x, i)
		for j in y:
			z = A[j, i]
			A[j, :] -= z*u/pivot
			I[j, :] -= z*v/pivot
	return I
