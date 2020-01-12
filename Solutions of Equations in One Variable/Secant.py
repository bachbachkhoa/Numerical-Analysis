# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from sympy import *

x = symbols('x')

def secant(f, a, b, err, Nmax):
	fa, fb = f.subs(x, a), f.subs(x, b)
	if fa*fb >= 0:
		return '\ninterval (a, b) wrong'
	else:
		if fa*diff(f, x, 2).subs(x, (a + b)/2) < 0:
			xBefore, d = a, b
		else:
			xBefore, d = b, a
		if abs(diff(f, x).subs(x, a)) < abs(diff(f, x).subs(x, b)):
			m, M = abs(diff(f, x).subs(x, a)), abs(diff(f, x).subs(x, b))
		else:
			M, m = abs(diff(f, x).subs(x, a)), abs(diff(f, x).subs(x, b))
		k = 0
		fd = f.subs(x, d)
		err = m/(M - m)*err
		while True:
			xAfter = xBefore - (d - xBefore)/(fd - f.subs(x, xBefore))*f.subs(x, xBefore)
			k += 1
			if abs(xAfter - xBefore) <= err or k >= Nmax:
				break
			else:
				xBefore = xAfter
	return float(xAfter), k
