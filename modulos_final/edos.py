import numpy as np
import matplotlib.pyplot as plt


def Euler(f,a,b,h,cond):
    n = int((b-a)/h)
    y_euler = [cond]
    for i in range(n):
        t_i = a+i*h
        y_euler.append(y_euler[i]+h*f(a+i*h,y_euler[i]))
    return np.linspace(a,b,n+1),y_euler

def Rk4(f,a,b,h,co):
    n = int((b-a)/h)
    w = [co]
    for i in range(n):
        k1 = h * f(a + i*h, w[i])
        k2 = h * f(a + i*h+h*0.5, w[i]+k1*0.5)
        k3 = h * f(a + i*h+h*0.5, w[i]+k2*0.5)
        k4 = h * f(a + i*h+h, w[i]+k3)
        w.append(w[i]+(1/6)*(k1+k2*2+2*k3+k4))
    return np.linspace(a,b,n+1),w