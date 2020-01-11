# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from sympy import *

x = symbols('x')

def Maclaurin(f, n):
	# f is symbolic, f is infinitely differentiable functions at 0
	result = []
	result.append(f.subs(x, 0))
	for i in range(1, n + 1):
		result.append(diff(f, x, i).subs(x, 0)/factorial(i))

	return result

def powerSeries(p, q, f, alpha, beta, n):
	# y'' + p(x)y' + q(x)y = f(x) and y(0) = alpha, y'(0) = beta
	# n is degree of y
	pM, qM, fM = Maclaurin(p, n), Maclaurin(q, n), Maclaurin(f, n)
	y = []
	y.append(alpha)
	y.append(beta)
	for i in range(2, n + 1):
		t = (fM[i - 2] - sum([(j + 1)*pM[i - j - 2]*y[j + 1] for j in range(i - 1)]) \
			- sum([qM[i - j - 2]*y[j] for j in range(i - 1)]))/i/(i - 1)
		y.append(t)

	return y
