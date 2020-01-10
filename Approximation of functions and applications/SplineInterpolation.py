# -*- coding: utf-8 -*-
from __future__ import unicode_literals

import numpy as np
from sympy import symbols, simplify

x = symbols('x')


## Tri Diagonal Matrix Algorithm(a.k.a Thomas algorithm) solver
def TDMAsolver(a, b, c, d):
	'''
	TDMA solver, a b c d can be NumPy array type or Python list type.
	refer to http://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
	and to http://www.cfd-online.com/Wiki/Tridiagonal_matrix_algorithm_-_TDMA_(Thomas_algorithm)
	'''
	nf = len(d) # number of equations
	ac, bc, cc, dc = map(np.array, (a, b, c, d)) # copy arrays
	for it in range(1, nf):
		mc = ac[it-1]/bc[it-1]
		bc[it] = bc[it] - mc*cc[it-1] 
		dc[it] = dc[it] - mc*dc[it-1]

	xc = bc
	xc[-1] = dc[-1]/bc[-1]

	for il in range(nf-2, -1, -1):
		xc[il] = (dc[il]-cc[il]*xc[il+1])/bc[il]

	return xc

def cubicSpline(data, k0, kn):
	# k0 = m0 = f''(x0), kn = mn = f''(xn)
	n = len(data) - 1
	X = data[:, 0]
	Y = data[:, 1]
	h = [X[i + 1] - X[i] for i in range(n)]
	d1 = [h[i] for i in range(1, n - 1)]
	d2 = [2*(h[i] + h[i + 1]) for i in range(n - 1)]
	d3 = [6*(h[i + 1]*Y[i] - (h[i + 1] + h[i])*Y[i + 1] + h[i]*Y[i + 2])/(h[i]*h[i + 1]) for i in range(n - 1)]
	d3[0] -= h[0]*k0
	d3[-1] -= h[-1]*kn
	m = TDMAsolver(d1, d2, d1, d3)
	m = m.tolist()
	m.insert(0, k0)
	m.insert(len(m), kn)

	a = [(m[i + 1] - m[i])/(6*h[i]) for i in range(n)]
	b = [m[i]/2 for i in range(n + 1)]
	c = [(Y[i + 1] - Y[i])/h[i] - (m[i + 1] + 2*m[i])*h[i]/6 for i in range(n)]
	S = [simplify(a[i]*(x - X[i])**3 + b[i]*(x - X[i])**2 + c[i]*(x - X[i]) + Y[i]) for i in range(n)]

	return S
