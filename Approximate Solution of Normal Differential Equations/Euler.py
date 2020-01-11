# -*- coding: utf-8 -*-
from __future__ import unicode_literals
# Author: Duynt

from sympy import *

x, y = symbols('x y')

def EulerForward(f, y0, a, b, h):
	n = int((b - a)/h)
	X = [a + i*h for i in range(0, n + 1)]
	result = zeros(n + 1, 2)
	result[:, 0] = X
	result[0, 1] = y0
	for i in range(1, n + 1):
		result[i, 1] = result[i - 1, 1] + h*f.subs(x, result[i - 1, 0]).subs(y, result[i - 1, 1])

	return result

def EulerBackward(f, y0, a, b, h):
	n = int((b - a)/h)
	X = [a + i*h for i in range(0, n + 1)]
	result = zeros(n + 1, 2)
	result[:, 0] = X
	result[0, 1] = y0
	for i in range(0, n):
		yt = result[i, 1] + h*f.subs(x, result[i, 0]).subs(y, result[i, 1])
		result[i + 1, 1] = result[i, 1] + 0.5*h*(f.subs(x, result[i, 0]).subs(y, result[i, 1]) \
							+ f.subs(x, result[i + 1, 0]).subs(y, yt))

	return result