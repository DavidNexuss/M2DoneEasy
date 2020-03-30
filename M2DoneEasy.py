from sympy import *

x = Symbol('x')
y = Symbol('y')
z = Symbol('z')
k = Symbol('k')
init_printing()

def taylor(function,x0,n,x = Symbol('x')):
    i = 0
    p = 0
    while i <= n:
        p = p + function.diff(x, i).subs(x, x0)/(factorial(i)) * (x - x0)**i
        i += 1
    return p
def lagrange(function,x0,n,x = Symbol('x')):
    return function.diff(x,n+1).subs(x,x0)/(factorial(n + 1)) * (x - x0)**(n + 1)
