# -*- coding: utf-8 -*-
from __future__ import unicode_literals
import numpy as np


'''Gauss-Jordan method to find root of linear system
   Input: augmented matrix of system
   Output: root of system
'''


def matRound(mat, decPts):
	'''Rounding elements of matrix'''
	for col in range(len(mat)):
		for row in range(len(mat[0])):
			mat[col][row] = round(mat[col][row], decPts)
	return mat

def convertMatrix(mat):
	'''Convert matrix to the form [[mat[0][0], 0, 0], [mat[0][1], 0, 1], ..., [mat[m][n], m, n]]'''
	result = []
	rowN = len(mat)
	colN = len(mat[0])
	for i in range(rowN):
		for j in range(colN):
			e = [mat[i][j], i, j]
			result.append(e)
	return result


def elementSolver(mat):

	# Select the element solver equal 1
	for i in mat:
		if abs(i[0]) == 1:
			return i

	# Select the element solver is integer and other than 0
	for i in mat:
		if isinstance(i[0], int) == True and abs(i[0]) != 1 and i[0] != 0:
			return i

	# Select the element solver is maximum element
	temp = []
	for i in mat:
		if i[0] != 0:
			j = i.copy()
			temp.append(j)
	maximum = temp[0]
	for i in temp:
		if abs(i[0]) > abs(maximum[0]):
			maximum = i
	return maximum



def gaussJordan(mat):
	'''Gauss-Jordan method to find root of linear system
	   Input: augmented matrix of system
	   Output: root of system
	'''

	A = [] # Coefficient matrix
	b = [] # Constant vector
	for i in mat:
		A.append(i[:-1])
		b.append(i[-1])

	rN = len(A)
	cN = len(A[0])

	convertA = convertMatrix(A)
	oldrIeS = []	# old idex row element solver
	oldcIeS = []	# old idex col element solver
	while len(convertA) != 0:
		eS = elementSolver(convertA)
		A_pq, p, q = eS
		oldrIeS.append(p)
		oldcIeS.append(q)
		for i in range(rN):
			if i != p:
				m = float(A[i][q])/float(A_pq)
				for j in range(cN):
					A[i][j] = A[i][j] - A[p][j]*m
				b[i] = b[i] - b[p]*m

		convertA = convertMatrix(A)
		convertA1 = []
		for i in convertA:
			if i[1] not in oldrIeS and i[2] not in oldcIeS:
				convertA1.append(i)
		convertA = convertA1


	A = matRound(A, 4)
	b = matRound(np.array(b).reshape(-1, 1).tolist(), 4)
	b = [item for elem in b for item in elem]

	# Convert coefficient matrix to form diagonal matrix
	for i in range(rN):
		for j in range(cN):
			if A[i][j] != 0 and i != j:
				A[i], A[j] = A[j], A[i]
				b[i], b[j] = b[j], b[i]


	# Find root of system
	x = []
	for i in range(rN):
		x.append(b[i]/A[i][i])

	return x



def main():
	inFile = input("\nEnter input file name: ")
	mat = np.loadtxt(inFile)
	mat = mat.tolist()

	# Find rank of coefficient matrix and augmented matrix
	n = len(mat[0]) - 1	# number of variables in the system
	A = []
	for i in mat:
		A.append(i[:-1])
	A = np.array(A)
	mat = np.array(mat)
	rmat = np.linalg.matrix_rank(mat)	# rank of augmented matrix
	rA = np.linalg.matrix_rank(A)		# rank of coefficient matrix

	# Solve system
	if rA == rmat == n:
		mat = mat.tolist()
		x = gaussJordan(mat)
		print("\nThe system has only a root:\n")
		print(x)
	elif rA == rmat < n:
		print("The system has infinite number of roots depending on", n - rA, "parameters")
	else:
		print("The system has no root")



if __name__ == '__main__':
	main()