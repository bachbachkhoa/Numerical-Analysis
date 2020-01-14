# Author: Duynt

import numpy as np

'''Gauss Elimination to find exact root of linear system Ax = b
   Input: augmented matrix of system
   Output: vector x or message: 'The system has no unique root'
'''

def gauss(A):
	'''Gauss Elimination to find exact root of linear system'''
	m = len(A)
	n = m + 1
	for i in range(m):
		# Find maximum pivot element
		pivots = [abs(A[k][i]) for k in range(i, m)]
		i_max = pivots.index(max(pivots)) + i
		if A[i_max][i] == 0:
			break
			return "The system has not unique root"
		else:
			# Convert matrix into upper triangular matrix
			A[i_max], A[i] = A[i], A[i_max]
			for k in range(i + 1, m):
				f = A[k][i]/A[i][i]
				for j in range(n):
					A[k][j] -= f*A[i][j]
				A[k][i] = 0
	# Solve equation Ax = b for an upper triangular matrix A    
	x = []
	for i in range(m - 1, -1, -1):
		x.insert(0, A[i][m]/A[i][i])
		for k in range(i - 1, -1, -1):
			A[k][m] -= A[k][i]*x[0]
	return x

def main():
	inFile = input("\nEnter input file name: ")
	A = np.loadtxt(inFile)
	A = A.tolist()

	# Find rank of coefficient matrix and augmented matrix
	n = len(A[0]) - 1	# number of variables in the system
	B = []
	for i in A:
		B.append(i[:-1])
	B = np.array(B)
	A = np.array(A)
	rA = np.linalg.matrix_rank(A)
	rB = np.linalg.matrix_rank(B)

	# Solve system
	if rA == rB == n:
		A = A.tolist()
		x = gauss(A)
		print("\nThe system has only a root:\n")
		print(x)
	elif rA == rB < n:
		print("The system has infinite number of roots depending on", n - rA, "parameters")
	else:
		print("The system has no root")


if __name__ == '__main__':
	main()
