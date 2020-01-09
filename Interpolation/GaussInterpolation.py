# -*- coding: utf-8 -*-
from __future__ import unicode_literals
import numpy as np

def differenceTable(data):
    # Create difference table
    n = len(data) - 1
    X = data[:, 0]
    Y = data[:, 1]

    f = np.zeros((n + 1, n + 1))
    f[:, 0] = Y
    for j in range(n):
        for i in range(n - j):
            f[i, j + 1] = f[i + 1, j] - f[i, j]

    return f

def findX0(x, data):
    # find start index
    X = data[:, 0]
    return np.where(abs(x - X) == min(abs(x - X)))[0][0]

def GaussForward(data, x):
    f = differenceTable(data)
    X = data[:, 0]
    startIndex = findX0(x, data)
    t = (x - X[startIndex])/(X[1] - X[0])
    deltaY, j = [], 0

    while True:
        if startIndex == 0 or startIndex + j > len(data) - 1:
            break
        deltaY.append(f[startIndex, j])
        if j%2 != 0:
            startIndex -= 1
        j += 1

    result = deltaY[0]
    deltaY = np.delete(deltaY, 0)
    i, prod = 2*[1]
    while i < len(deltaY) + 1:
        if i%2 != 0:
            prod *= (t + i//2)
        else:
            prod *= (t - i//2)
        result += prod*deltaY[i - 1]/np.math.factorial(i)
        i += 1

    return result

def GaussBackward(data, x):
    f = differenceTable(data)
    X = data[:, 0]
    startIndex = findX0(x, data)
    t = (x - X[startIndex])/(X[1] - X[0])
    deltaY, j = [], 0

    while True:
        if (startIndex == 0 and j%2 == 0) or (startIndex + j == len(data) - 1 and j%2 != 0):
            break
        deltaY.append(f[startIndex, j])
        if j%2 == 0:
            startIndex -= 1
        j += 1

    result = deltaY[0]
    deltaY = np.delete(deltaY, 0)
    i, prod = 2*[1]
    while i < len(deltaY) + 1:
        if i%2 == 0:
            prod *= (t + i//2)
        else:
            prod *= (t - i//2)
        result += prod*deltaY[i - 1]/np.math.factorial(i)
        i += 1

    return result

def Stirling(data, x):
    f = differenceTable(data)
    X = data[:, 0]
    startIndex = findX0(x, data)
    t = (x - X[startIndex])/(X[1] - X[0])
    deltaY, j = [], 0

    while True:
        if (startIndex == 0 and j%2 != 0) or startIndex + j > len(data) - 1:
            break
        if j%2 == 0:
            deltaY.append(f[startIndex, j])

        if j%2 != 0:
            deltaY.append([f[startIndex, j], f[startIndex - 1, j]])
            startIndex -= 1
        j += 1

    result = deltaY[0] + t*(deltaY[1][0] + deltaY[1][1])/2
    deltaY = np.delete(deltaY, [0, 1])
    i, prod1, prod2 = 2, t, 1
    while i < len(deltaY) + 2:
        if i%2 == 0:
            prod2 *= (t**2 - (i//2 - 1)**2)
            result += prod2*deltaY[i - 2]/np.math.factorial(i)
        else:
            prod1 *= (t**2 - (i//2)**2)
            result += prod1*(deltaY[i - 2][0] + deltaY[i - 2][1])/(2*np.math.factorial(i))
        i += 1

    return result

def Bessel(data, x):
    f = differenceTable(data)
    X = data[:, 0]
    startIndex = findX0(x, data)
    t = (x - X[startIndex])/(X[1] - X[0])
    deltaY, j = [], 0

    while True:
        if startIndex < 0 or (startIndex + j == len(data) - 1 and j%2 == 0):
            break
        if j%2 != 0:
            deltaY.append(f[startIndex, j])
            startIndex -= 1
        if j%2 == 0:
            deltaY.append([f[startIndex, j], f[startIndex + 1, j]])
        j += 1

    result = (deltaY[0][0] + deltaY[0][1])/2 + (t - 1/2)*deltaY[1]
    deltaY = np.delete(deltaY, [0, 1])
    i, prod1, prod2 = 2, t - 1/2, 1
    while i < len(deltaY) + 2:
        if i%2 == 0:
            prod2 *= (t + i//2 - 1)*(t - i//2)
            result += prod2*(deltaY[i - 2][0] + deltaY[i - 2][1])/(2*np.math.factorial(i))
        else:
            prod1 *= (t + i//2 - 1)*(t - i//2)
            result += prod1*deltaY[i - 2]/np.math.factorial(i)
        i += 1

    return result

