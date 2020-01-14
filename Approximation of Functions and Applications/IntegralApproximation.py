
# from sympy import *
# from scipy import optimize
# import numpy as np

# x = symbols('x')

def trapezoidalRule(f, a, b, n):
	# f2 = diff(f, x, 2)
	# y = lambdify(x, -abs(f2))
	# f = lambdify(x, f)
	# xmaxf2 = optimize.minimize_scalar(y, bounds=[a, b], method='bounded')['x']
	# maxf2 = abs(f2.subs(x, xmaxf2))
	# n = int(sqrt(maxf2*(b - a)**3/12/err)) + 1
	h = (b - a)/n
	return h/2*(f(a)+ f(b) + 2*sum([f(a + i*h) for i in range(1, n)]))

def SimpsonRule(f, a, b, n):
	if n%2 != 0:
		return 'n must be an even number'
	else:
		h = (b - a)/n
		return h/3*(f(a)+ f(b) + 4*sum([f(a + i*h) for i in range(1, n, 2)]) + 2*sum([f(a + i*h) for i in range(2, n - 1, 2)]))
