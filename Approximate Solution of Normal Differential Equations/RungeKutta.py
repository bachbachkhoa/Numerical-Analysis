# -*- coding: utf-8 -*-
from __future__ import unicode_literals
# Author: Duynt

from sympy import *

x, y = symbols('x y')

def RungeKutta3(f, y0, a, b, h):
	n = int((b - a)/h)
	result = zeros(n + 1, 2)
	result[:, 0] = [a + i*h for i in range(0, n + 1)]
	result[0, 1] = y0
	for i in range(n):
		k1 = h*f.subs(x, result[i, 0]).subs(y, result[i, 1])
		k2 = h*f.subs(x, result[i, 0] + h/2).subs(y, result[i, 1] + k1/2)
		k3 = h*f.subs(x, result[i, 0] + h).subs(y, result[i, 1] - k1 + 2*k2)
		result[i + 1, 1] = result[i, 1] + (k1 + 4*k2 + k3)/6

	return result

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
