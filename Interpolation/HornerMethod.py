# -*- coding: utf-8 -*-
from __future__ import unicode_literals


'''Horner's method and its applications '''

# Horner's method
def Horner(a, x0):
	'''Input: Polynomial a and x0 (a is a list of coefficients)
	   Output: a(x0) and quotient of the division polynomial a for (x - x0)
	'''
	n = len(a) - 1			# degree of the polynomial is n
	b = []					# b\{b[n]} is a polynomial satisfying: a = (x - x0)*b + a(x0)
	b.insert(len(b), a[0])	# b[n] is value of the polynomial at x0

	for i in range(1, n + 1):
		t = a[i] + b[i - 1]*x0
		b.insert(len(b), t)

	a_x0 = b[-1]
	b.pop()
	return a_x0, b
  
# Applications of Horner's method

