# -*- coding: utf-8 -*-
from __future__ import unicode_literals
import numpy as np

A = [[1, 1, -3, 2], [1, -2, 0, -1], [0, 1, 1, 3], [2, -3, 2, 0]]
b = [6, -6, 16, 6]

def isIdentity(mat):
	for row in range(len(mat)):
		for col in range(len(mat[0])):
			if row == col and mat[row][col] != 1:
				return False 
			elif row != col and mat[row][col] != 0:
				return False
	return True

def convertMatrix(mat):
	result = []
	rowN = len(mat)
	colN = len(mat[0])
	for i in range(rowN):
		for j in range(colN):
			e = [mat[i][j], i, j]
			result.append(e)
	return result

def elementSolver(cmat):

	# Select the element solver equal 1
	for i in cmat:
		if abs(i[0]) == 1:
			return i

	# Select the element solver is integer and other than 0
	for i in cmat:
		if isinstance(i[0], int) == True and abs(i[0]) != 1 and i[0] != 0:
			return i

	# Select the element solver is maximum element
	temp = []
	for i in cmat:
		if i[0] != 0:
			j = i.copy()
			temp.append(j)
	maximum = temp[0]
	for i in temp:
		if abs(i[0]) > abs(maximum[0]):
			maximum = i
	return maximum



def gaussJordan(A, b):
	if len(A) != len(b):
		return 'A and b are not conformable'
	else:
		convertA = convertMatrix(A)
		oldrIeS = []	# old idex row element solver
		oldcIeS = []	# old idex col element solver
		while len(convertA) != 0:
			eS = elementSolver(convertA)
			A_pq, p, q = eS
			oldrIeS.append(p)
			oldcIeS.append(q)
			for i in range(0, len(A)):
				if i != p:
					m = A[i][q]/A_pq
					for j in range(0, len(A[0])):
						A[i][j] = A[i][j] - A[p][j]*m
					b[i] = b[i] - b[p]*m

			convertA = convertMatrix(A)
			convertA1 = []
			for i in convertA:
				if i[1] not in oldrIeS and i[2] not in oldcIeS:
					convertA1.append(i)
			convertA = convertA1

		x = []
		for i in range(len(A)):
			for j in range(len(A[0])):
				if A[i][j] != 0:
					x.append(b[i]/A[i][j])
		return x

print(gaussJordan(A, b))