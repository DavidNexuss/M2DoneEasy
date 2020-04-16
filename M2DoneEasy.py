from sympy import *

x = Symbol('x')
y = Symbol('y')
z = Symbol('z')
k = Symbol('k')
c = Symbol('c')
init_printing()

#Per utsar taylor es recomanable que la funció no utilitze fraccionaris
#de la forma 231/243 en lloc, Rational(231,243).
#Sinó simpy evaluarà el resultat en decimal.
#Es pot usar expand per expandir la expressió donada per la funció i així que 
#sympy la imprimeixi amb el format convencional
#dels polinomis de taylor

def taylor(function,x0,n,x = Symbol('x')):
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
print("Benvingut/da a la biblioteca de funcions de M2 de +K linux UPC")
