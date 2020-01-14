'''
   *****************************************
   * Classic Iterative Method to find      *
   * inverse of diagonally dominant matrix *
   *****************************************
   *     Author     *        Duynt         *
   *****************************************
'''

'''
Usage: classicIteration.py -i <inputfile> -o <outputfile> -x <initialfile> -e <error>
Help:  classicIteration.py -h
'''


import numpy as np
import sys, getopt

def isSquare(m):
	'''Check square matrix. Return true if the matrix is ​​a square matrix '''
	return all (len (row) == len (m) for row in m)

def isDDM(A):
	'''Check diagonally dominant matrix. Return
	   0: Not diagonally dominant matrix
	   1: Row diagonal dominance
	   2: Column diagonal dominance'''
	d = len(A)
	rowCheck = 1
	colCheck = 2
	for i in range(d):
		rowSum = 0
		colSum = 0
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

def diagonalMatrixT(A):
	'''Create diagonal matrix T and calculate T^-1'''
	d = len(A)
	T = np.eye(d)
	invT = np.eye(d)
	for i in range(d):
		T[i][i] = 1/A[i][i]
	for i in range(d):
		invT[i][i] = A[i][i]
	return (T, invT)

def rowDDM(A, T, X, err):
	'''Row diagonal dominance'''
	E = np.eye(len(A))
	X1 = np.dot(np.subtract(E, np.dot(T, A)), X) + np.dot(T, E)
	k = 0
	norm = np.linalg.norm(np.subtract(E, np.dot(T, A)), np.Inf)
	err = err*(1 - norm)/norm
	m = np.linalg.norm(np.subtract(X1, X), np.Inf)
	X_before = X1
	while m > err:
		k += 1
		X_after = np.dot(np.subtract(E, np.dot(T, A)), X_before) + np.dot(T, E)
		m = np.linalg.norm(np.subtract(X_after, X_before), np.Inf)
		X_before = X_after
	return (X_before, k)

def colDDM(A, T, Y, err):
	'''Column diagonal dominance'''
	E = np.eye(len(A))
	Y1 = np.dot(np.subtract(E, np.dot(A, T)), Y) + E
	k = 0
	norm = np.linalg.norm(np.subtract(E, np.dot(A, T)), 1)
	err = err*(1 - norm)/(norm*np.linalg.norm(T, 1))
	m = np.linalg.norm(np.subtract(Y1, Y), 1)
	Y_before = Y1
	while m > err:
		k += 1
		Y_after = np.dot(np.subtract(E, np.dot(A, T)), Y_before) + E
		m = np.linalg.norm(np.subtract(Y_after, Y_before), 1)
		Y_before = Y_after
	return (np.dot(T, Y_before), k)


def main(argv):
	try:
		opts, args = getopt.getopt(argv,"hi:o:x:e:",["ifile=","ofile=","xfile=","err="])
	except getopt.GetoptError:
		print('classicIteration.py -i <inputfile> -o <outputfile> -x <initialfile> -e <error>')
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print('-h: help')
			print('-i: input file name')
			print('-o: output file name')
			print('-x: initial file name')
			print('-e: error')
			sys.exit(0)
		elif opt in ("-i", "--ifile"):
			inFile = arg
		elif opt in ("-o", "--ofile"):
			outFile = arg
		elif opt in ("-x", "--xfile"):
			initialFile = arg
		elif opt in ("-e", "--err"):
			err = float(arg)
	try:
		A = np.loadtxt(inFile)
		X = np.loadtxt(initialFile)
		checkDDM = isDDM(A)
		if isSquare(A) == False:
			print("\nThe matrix is not square matrix")
		elif len(X) != len(A) or len(X[0]) != len(A):
			print("\nThe initialization matrix and the input matrix are not the same size")
		elif np.linalg.det(A) == 0:
			print("\nThe matrix is not nonsingular matrix")
		elif checkDDM == 0:
			print("\nThe matrix is not diagonally dominant matrix")
		else:
			T, invT = diagonalMatrixT(A)
			if checkDDM == 1:
				result1, nIteration1 = rowDDM(A, T, X, err)
				NoI1 = "\n\nNumber of iteration: " + str(nIteration1)
				np.savetxt(outFile ,result1)
				with open(outFile, 'a') as f:
					f.write(NoI1)
			elif checkDDM == 2:
				Y = np.dot(invT, X)
				result2, nIteration2 = colDDM(A, T, Y, err)
				NoI2 = "\n\nNumber of iteration: " + str(nIteration2)
				np.savetxt(outFile ,result2)
				with open(outFile, 'a') as f:
					f.write(NoI2)
	except:
		print("An error occurred during processing")


if __name__ == '__main__':
	main(sys.argv[1:])
