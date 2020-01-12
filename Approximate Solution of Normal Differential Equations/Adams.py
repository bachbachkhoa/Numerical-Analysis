# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from sympy import *
# import numpy as np

x, y = symbols('x y')

def RungeKutta4(f, y0, a, b, h):
	n = int((b - a)/h)
	result = zeros(n + 1, 2)
	result[:, 0] = [a + i*h for i in range(0, n + 1)]
	result[0, 1] = y0
	for i in range(n):
		k1 = h*f.subs(x, result[i, 0]).subs(y, result[i, 1])
		k2 = h*f.subs(x, result[i, 0] + h/2).subs(y, result[i, 1] + k1/2)
		k3 = h*f.subs(x, result[i, 0] + h/2).subs(y, result[i, 1] + k2/2)
		k4 = h*f.subs(x, result[i, 0] + h).subs(y, result[i, 1] + k3)
		result[i + 1, 1] = result[i, 1] + (k1 + 2*k2 + 2*k3 + k4)/6

	return result

def AdamsInterpolation(f, y0, a, b, h, err):
	n = int((b - a)/h)
	result = zeros(n + 1, 2)
	result[:, 0] = [a + i*h for i in range(0, n + 1)]
	rk4 = RungeKutta4(f, y0, a, b, h)
	result[:3, 1] = rk4[:3, 1]
	for i in range(2, n):
		dyi = f.subs(x, result[i, 0]).subs(y, result[i, 1])
		dyi1 = f.subs(x, result[i - 1, 0]).subs(y, result[i - 1, 1])
		dyi2 = f.subs(x, result[i - 2, 0]).subs(y, result[i - 2, 1])
		sigma = h/24*(19*dyi - 5*dyi1 + dyi2)
		yBefore = result[i, 1] + 9/24*h*f.subs(x, result[i + 1, 0]).subs(y, result[i, 1]) + sigma
		while True:
			yAfter = result[i, 1] + 9/24*h*f.subs(x, result[i + 1, 0]).subs(y, yBefore) + sigma
			if abs(yAfter - yBefore) <= err:
				break
			else:
				yBefore = yAfter
		result[i + 1, 1] = yAfter

	return result

# ad = np.array(ad.tolist()).astype(np.float64)
# np.savetxt('a.txt', ad)
