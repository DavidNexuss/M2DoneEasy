from sympy import *

x = Symbol('x')
y = Symbol('y')
z = Symbol('z')
k = Symbol('k')
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
def lagrange(function,x0,n,x = Symbol('x')):
    return function.diff(x,n+1).subs(x,x0)/(factorial(n + 1)) * (x - x0)**(n + 1)


print("Benvingut/da a la biblioteca de funcions de M2 de +K linux UPC")
