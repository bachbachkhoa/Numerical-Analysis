# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from sympy import *
from scipy import optimize

'''Solution of equation: x = phi(x)
   phi(x) in [a, b], abs(phi'(x)) < 1 with all x in [a, b]'''

x = symbols('x')

def fixedPoint(phi, a, b, x0, err, Nmax):
	# phi is symbolic
	f = x - phi
	if f.subs(x, a)*f.subs(x, b) >= 0:
		return '\ninterval (a, b) wrong'
	else:
		if x0 < a or x0 > b:
			return '\nx0 must be in [{0}, {1}]'.format(a, b)
		else:
			phi1 = diff(phi, x)
			y = lambdify(x, -abs(phi1))
			xmax = optimize.minimize_scalar(y, bounds=[a, b], method='bounded')['x']
			q = phi1.subs(x, xmax)
			err *= (1 - q)/q
			k = 0
			while True:
				x1 = phi.subs(x, x0)
				k += 1
				if abs(x1 - x0) <= err or k >= Nmax:
					break
				else:
					x0 = x1
		return float(x1), k
