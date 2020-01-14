# Author: Duynt

import numpy as np


'''Cholesky method to find root of linear system
   Input: augmented matrix of system (coefficient matrix is symmetric matrix)
   Output: root of system
'''


def Cholesky(A, b):
	
	# A is symmetric matrix
	# Decomposite A = L.T*L such that L is lower triangular matrix
	L = np.zeros_like(A)
	n = len(L)
	for i in range(n):
		for j in range(i+1):
			if i==j:
				val = A[i][i] - np.sum(np.square(L[i][:i]))
				L[i][i] = np.sqrt(val)
			else:
				L[i][j] = (A[i][j] - np.sum(L[i][:j]*L[j][:j]))/L[j][j]

	# Calculate y such that L*y = b
	y = []
	y.insert(0, b[0]/L[0][0])
	for i in range(1, n):
		b[i] -= np.sum(L[i][:i]*y[:i])
		y.insert(len(y), b[i]/L[i][i])

	# Calculate x such that L.T*x = y
	x = []
	x.insert(0, y[-1]/L[-1][-1])
	for i in range(n - 2, -1, -1):
		for k in range(i + 1, n):
			y[i] -= L[k][i]*x[k - n]
		x.insert(0, y[i]/L[i][i])
	return x

def main():
	inFile = input("\nEnter input file name: ")
	m = np.loadtxt(inFile)	# Augmented matrix
	m = m.tolist()

	A = []	# Coefficient matrix (symmetrix matrix)
	b = []	# Constant vector
	for i in m:
		A.append(i[:-1])
		b.append(i[-1])
	rm = np.linalg.matrix_rank(np.array(m))	# Rank of augmented matrix
	rA = np.linalg.matrix_rank(np.array(A))	# Rank of coefficient matrix

	# Solve system
	if rA == rm == len(b):
		A = np.array(A).astype(complex)
		b = np.array(b).astype(complex)
		A1 = np.dot(A.T, A)
		b = np.dot(A.T, b)
		x = Cholesky(A1, b)
		print("\nThe system has only a root:\n")
		print(x)
	
	elif rA == rm < len(b):
		print("The system has infinite number of roots depending on", len(b) - rA, "parameters")
	
	else:
		print("The system has no root")



if __name__ == '__main__':
	main()
