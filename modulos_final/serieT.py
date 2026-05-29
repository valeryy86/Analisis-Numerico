import sympy as sp
import numpy as np
from math import factorial

x = sp.symbols('x')

def taylor_serie(f,x0,n):
    polinomio = 0
    for k in range(n+1):
        df = sp.diff(f,x,k)
        df_eval = sp.lambdify(x,df)
        polinomio += (df_eval(x0)/factorial(k))*(x - x0)**k
    return sp.expand(polinomio)

def cota_truncamiento(f,n,x0,x_point):
    df = sp.diff(f,x,n+1)
    df_eval = sp.lambdify(x,df)
    x_array = np.linspace(min(x0,x_point), max(x0,x_point),200)
    maximo = max(abs(df_eval(x_array)))
    cota = maximo*(abs(x0-x_point))**(n+1)/factorial(n+1)
    return cota