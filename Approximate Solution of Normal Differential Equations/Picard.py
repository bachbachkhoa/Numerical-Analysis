# Author: Duynt

from sympy import *

x, y = symbols('x y')


def Picard(f, x0, y0, n):
	'''f is symbolic, y(x0) = y0, n is number of iterations'''
	k = 0
	yBefore = y0
	while True:
		yAfter = y0 + integrate(f.subs(y, yBefore), (x, x0, x))
		k += 1
		if k >= n:
			break
		else:
			yBefore = yAfter

	return yAfter
