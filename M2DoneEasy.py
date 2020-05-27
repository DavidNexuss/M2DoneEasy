from sympy import *
from sympy.vector import CoordSys3D, Del

x0 = Symbol('x0')
x1 = Symbol('x1')
x2 = Symbol('x2')
x3 = Symbol('x3')
x4 = Symbol('x4')

x = Symbol('x')
y = Symbol('y')
z = Symbol('z')
k = Symbol('k')
c = Symbol('c')
R = CoordSys3D('R')

a = Symbol('a')
b = Symbol('b')
c = Symbol('c')
d = Symbol('d')
e = Symbol('e')
f = Symbol('f')

auxiliars = [
Symbol('a'),
Symbol('b'),
Symbol('c'),
Symbol('d'),
Symbol('e'),
Symbol('f')
]

n = Symbol('n')
delop = Del()
init_printing()

#Per utsar taylor es recomanable que la funció no utilitze fraccionaris
#de la forma 231/243 en lloc, Rational(231,243).
#Sinó simpy evaluarà el resultat en decimal.
#Es pot usar expand per expandir la expressió donada per la funció i així que 
#sympy la imprimeixi amb el format convencional
#dels polinomis de taylor

#def mul_subst(function,vector,variables):
#    for x in vector:
#        function = function.subs(x,
#def taylor_mul(function,a,k):

def taylor_mul_1(f,a,b):
    return f.subs(x,a).subs(y,b) + f.diff(x).subs(x,a).subs(y,b) * (x - a) + f.diff(y).subs(x,a).subs(y,b) * (y - b)
def taylor_mul_2(f,a,b):
    return taylor_mul_1(f,a,b) + (1/2)*(f.diff(x,2).subs(x,a).subs(y,b)*(x-a)**2 + 2*f.diff(x).diff(y).subs(x,a).subs(y,b)*(x-a)*(y-b)+ f.diff(y,2).subs(x,a).subs(y,b)*(y-b)**2)
def mul_vector(a,b):
    c = 0
    for i in range(len(a)):
        c += a[i] * b[i]
    return c
def restrict(F,restriction):
    for i in range(len(restriction)):
        F = F + auxiliars[i]*restriction[i]
    return F
def intersection_matrix(v,w):
    v = v.row_join(v)
    w = w.row_join(w * 0)
    v = v.col_join(w)
    return v
#reuires the intersection matrix m of vector spaces v and w and computes its vector space intersection
def Zassenhaua(m):
    cols = int(m.cols / 2)
    rref = m.rref()
    r = 0
    for i in rref[1]:
        if i < cols:
            r = r + 1
    for i in range(r):
        rref[0].row_del(0)
    for i in range(cols):
        rref[0].col_del(0)
    return rref[0]
def Jacobian(functions,variables):
    mat = []
    for f in functions:
        row = []
        for v in variables:
            row.append(f.diff(v))
        mat.append(row)
    return Matrix(mat)
def Hessian(function,variables):

    mat = []
    for i in variables:
        row = []
        for j in variables:
            row.append(function.diff(i).diff(j))
        mat.append(row)
    return Matrix(mat)
def Gradient(p_function, variables = 0):
    row = []
    if variables == 0:
        variables = p_function.free_symbols
    for i in variables:
        row.append(p_function.diff(i))
    return Matrix(row)
def generic_taylor(f):
    return Sum(f.diff((x,n)).subs(x,x0)/factorial(n) * (x - x0)**n,(n,0,n))

def taylor(function,_x0,_n,x = Symbol('x')): 
    i = 0
    p = 0
    while i <= n:
        p = p + function.diff(x, i).subs(x, x0)/(factorial(i)) * (x - x0)**i
        i += 1
    return p

#Retorna l'expressió de l'error en forma de lagrange
def lagrange(function,x0,xf,n,x = Symbol('x')):
    return function.diff(x,n+1).subs(x,c)/(factorial(n + 1)) * (UnevaluatedExpr(xf) - UnevaluatedExpr(x0))**(n + 1)

def bi(f,a,b,it,x= Symbol('x')):
    numbers = []
    numbers.append(a)
    while it != 0:
        c = (a + b) / 2
        numbers.append(c)
        A = f.subs(x,a)
        B = f.subs(x,b)
        C = f.subs(x,c)
        if C == 0: break
        if C * A > 0:
            a = c
        else:
            b = c
        it = it - 1;
    return numbers
#Retorna les iteracions necessaries del metode de la biseccio per a una precisio desitjada
def bi_precision(b,a,c):
    return ceiling(ln((a - b)/c)/ln(2))
def secante(f,x0,x1,it,x = Symbol('x')):
    numbers = []
    numbers.append(x0)
    numbers.append(x1)
    while it != 0:
        n = len(numbers)
        xn = numbers[n-1]
        xn1 = numbers[n-2]
        dis = f.subs(x,xn) - f.subs(x,xn1)
        if dis != 0:
            numbers.append(N(
                xn - ((xn - xn1)/(f.subs(x,xn) - f.subs(x,xn1))) * f.subs(x,xn)
            ))
        else:
            break
        it = it - 1
    return numbers
#Retorna el valor aproximat d'igualar a 0 la funció f amb el métode newton rhapson
def nr(f,x0,it,x = Symbol('x')):
    g = diff(f,x)
    numbers = []
    numbers.append(x0)
    while it != 0:
        x0 = x0 - N(f.subs(x,x0)/g.subs(x,x0))
        numbers.append(x0)
        it = it - 1
    return numbers

#Retorna el valor aproximat d'igualar a 0 la funció f amb el métode newton rhapson en funcio de una precisió p
def nr_p(f,x0,p,x = Symbol('x')):
    g = diff(f,x)
    numbers = []
    x1 = x0
    numbers.append(x1)
    while Abs(x1 - x0) >= p or Abs(f.subs(x,x1)) >= p:
        t = x1
        x1 = x0 - N(f.subs(x,x0)/g.subs(x,x0))
        numbers.append(x1)
        x0 = t
    return numbers
def suc(f,x0,it,x = Symbol('x')):
    numbers = []
    numbers.append(x0)
    while it != 0:
        x0 = f.subs(x,x0)
        numbers.append(x0)
        it = it - 1
    return numbers

#-----------------------------------------Algebra Lineal-------------------------------


print("Benvingut/da a la biblioteca de funcions de M2 y M1 de +K linux UPC")
