'''
Newton method to find root of equation 
Author: Duynt
'''

from sympy import *
import numpy as np

'''
Condition: 
* f, f', f'' are continuous functions defined on the interval [a, b], with f (a) and f (b)
of opposite sign. Root of f in this interval is unique.
* f' have the same sign for every x in (a, b) and f'' too.
'''

x = Symbol('x')

def Newton(f, a, b, err, Nmax):
    if f.subs(x, a)*f.subs(x, b) >= 0:
        return '\ninterval (a, b) wrong'
    else:
        df = f.diff(x)
        m = min(abs(df.subs(x, a)), abs(df.subs(x, b)))
        if f.subs(x, a)*f.diff(x, 2).subs(x, a) > 0:
            x0 = a
        else:
            x0 = b
        k = 0
        while abs(f.subs(x, x0))/m > err and k < Nmax:
            x0 = x0 - f.subs(x, x0)/df.subs(x, x0)
            k = k + 1
        x0 = float(x0)
        return (x0, k)


def main():
    f = x**3 + x - 1000
    if Newton(f, 9, 10, 1e-5, 100) == '\ninterval (a, b) wrong':
        print('\ninterval (a, b) wrong')
    else:
        root, k = Newton(f, 9, 10, 1e-5, 100)
        print("\nroot of equation: ", root)
        print("\nnumber of iteration: ", k)

if __name__ == '__main__':
    main()
    