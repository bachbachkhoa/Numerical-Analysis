import numpy as np

#priori evaluation
def bisection_e(f, a, b, err):
	if f(a) * f(b) >= 0:
		print("Bisection method fail!")
		return None
	a_n = a
	b_n = b
	N = int(np.log((b_n - a_n)/err)/np.log(2) - 1) + 1
	for n in range(1, N + 1):
		m_n = (a_n + b_n)/2
		f_m_n = f(m_n)
		if f(a_n)*f_m_n < 0:
			b_n = m_n
		elif f(b_n)*f_m_n < 0:
			a_n = m_n
		elif f_m_n == 0:
			print("Found root of f!")
			return m_n
		else:
			print("Bisection method fail!")
			return None
	return (a_n + b_n)/2

#posteriori evaluation
def bisection_l(f, a, b, err):
	if f(a)*f(b) > 0:
		print("Bisection method fail!")
		return None
	a_0 = a
	b_0 = b
	while np.abs(b - a) > err:
		c = (a_0 + b_0)/2
		z = f(c)
		if z*f(a) < 0:
			b_0 = c
		elif z*f(b) < 0:
			a_0 = c
		elif z == 0:
			print("Found root of f!")
			return c
		else:
			print("Bisection method fail!")
			return None
	return (a_0 + b_0)/2



f = lambda x: x*x - x -1
root_e = bisection_e(f, 1, 2, pow(10, -9))
root_l = bisection_l(f, 0, 2, pow(10, -9))
print(root_e)
print(root_l)
