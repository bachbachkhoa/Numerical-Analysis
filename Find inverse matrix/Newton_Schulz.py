'''
   **************************************************
   * Newton Schulz Method to find inverse of matrix *
   **************************************************
   *         Author        *          Duynt         *
   **************************************************
'''


import numpy as np

def isSquare(A):
    '''Check square matrix'''
    return all (len (row) == len (A) for row in A)

'''Initialize X0'''
def initialX(A):
    alpha = np.linalg.norm(A, 2)**-2
    return np.dot(alpha, np.transpose(A))

def NewtonSchulz(A, X, err):
    '''Newton Schulz method'''
    d = len(A)
    E = np.eye(d)
    G0 = np.subtract(E ,np.dot(A, X))
    normG0 = np.linalg.norm(G0, 2)
    normInitialX = np.linalg.norm(X, 2)
    k = 0
    m = normG0**(2**k)
    err = err*(1 - normG0)/normInitialX
    while m > err:
        X = np.dot(X, np.subtract(np.dot(2,E), np.dot(A, X)))
        k = k + 1
        m = normG0**(2**k)
    return (X, k)

def main():
    try:
        inFile = input("\nEnter input file name: ")
        A = np.loadtxt(inFile)
        if isSquare(A) == False:
            print("\nThe matrix is not a square matrix")
        elif np.linalg.det(A) == 0:
            print("\nThe matrix is not an invertible matrix")
        else:
            outFile = input("\nEnter output file name: ")
            err = float(input("\nEnter error: "))
            YorN = input('\nEnter X0 from file?(y/n) ')
            while YorN != 'y' and YorN != 'n' and YorN != 'Y' and YorN != 'N':
                YorN = input('\ny or n??? Retype: ')
            if YorN == 'y' or YorN == 'Y':
                initialFile = input("\nEnter initial file name: ")
                X = np.loadtxt(initialFile)
                if len(X) != len(A) or len(X[0]) != len(A):
                    print("\nThe initialization matrix and the input matrix are not the same size")
                else:
                    invA, k = NewtonSchulz(A, X, err)
                    np.savetxt(outFile, invA)
                    with open(outFile, 'a') as f:
                        f.write("\n\nNumber of iteration: " + str(k))
            elif YorN == 'n' or YorN == 'N':
                X = initialX(A)
                invA, k = NewtonSchulz(A, X, err)
                np.savetxt(outFile, invA)
                with open(outFile, 'a') as f:
                    f.write("\n\nNumber of iteration: " + str(k))
    except:
        print('An error occurred during processing')

if __name__ == '__main__':
	main()
