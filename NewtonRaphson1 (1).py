import sympy as sp
from sympy.utilities.lambdify import lambdify
import termcolor


def f(x):
    """"" Polynomial Function Definietion"""
    return x **2 +7*x-5


def g(y):
    """"" Defining Derviative Of A Function"""
    x = sp.symbols('x')
    f = x ** 2 + 7 * x - 5
    fp = f.diff(x)
    fp = lambdify(x, fp)
    return fp(y)


def NewtonRaphson(xs, e):
    """"" Defining The Newthon Raphson Method."""""
    print(termcolor.colored('\n\n~~Newton Raphson Method~~','red'))
    counter = 1
    term = True
    while term:
        if g(xs) == 0.0:
            print('Cannot Divide By Zero!')
            break
        else:
            xi = xs - f(xs) / g(xs)
            print('Iteration: %d, x1 = %0.6f , f(x1) = %0.6f' % (counter, xi, f(xi)))
            xs = xi
            counter+= 1
            term = abs(f(xi)) > e

    print('Root is: %0.8f' % xi)

x1 = int(input(' Please Enter Range: '))
x2 = int(input('Please Enter Range: '))
x0 = (x2 - x1) / 2

# Convert To Float
x0 = float(x0)
e = float(0.00001)
# Activate Newton Raphson Method:
NewtonRaphson(x0, e)