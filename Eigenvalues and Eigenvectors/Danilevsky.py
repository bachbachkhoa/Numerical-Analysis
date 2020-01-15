# -*- coding: utf-8 -*-
from __future__ import unicode_literals
# Author: Duynt

import numpy as np

def checkFrobeniusMatrix(a):
	n = len(a)
	x = a[1:, :n - 1].copy()
	y = a[1:, n - 1].copy()
	if np.all(np.equal(x, np.eye(n - 1))) and all(y == 0):
		return True
	return False

def Danilevski(a):
	if checkFrobeniusMatrix(a):
		p = a[0, :].copy()
		p = (-1)**(n + 1)*np.insert(p, 0, -1)
		return p
	else:
		n = len(a)
		result = 1
		for k in range(n - 1):
			if a[n - k - 1, n - k - 2] != 0:
				m, invm = np.eye(n), np.eye(n)
				m[n - k - 2, :] = a[n - k - 1, :].copy()
				invm[n - k - 2, :] =  -a[n - k - 1, :]/a[n - k - 1, n - k - 2]
				invm[n - k - 2, n - k - 2] = 1/a[n - k - 1, n - k - 2]
				a = np.dot(np.dot(m, a), invm)
			elif all(a[n - k - 1, :n - k - 1] == 0) == True:
				a3 = a[n - k - 1:, n - k - 1:].copy()
				a1 = a[:n - k - 1, :n - k - 1].copy()
				p3 = a3[0, :].copy()
				p3 = (-1)**(len(a3) + 1)*np.insert(p3, 0, -1)
				return np.polymul(np.polymul(result, p3), Danilevski(a1))
			elif any(a[n - k - 1, :n - k - 1] != 0) == True:
				swapPos = np.where(a[n - k - 1, :n - k - 1] != 0)[0][0]
				a[:, [swapPos, n - k - 2]] = a[:, [n - k - 2, swapPos]]
				a[[swapPos, n - k - 2], :] = a[[n - k - 2, swapPos], :]
				return Danilevski(a)
		
		p1 = a[0, :].copy()
		p1 = (-1)**(n + 1)*np.insert(p1, 0, -1)
		result = np.polymul(result, p1)

		return result

a = np.loadtxt('a.txt')
print(Danilevski(a))