#!/usr/bin/env python3
from scipy import *
from scipy import integrate
import numpy as np

omega = 0.0228
E = 0.0655
ionpot = 0.579

n = 36.0

ac = 0.001
q = 1.0

t_start = -pi/(2*omega)
t_end = 3*pi/(2*omega)
v0_min = 1.0
v0_max = -1.0

def El(t):
    return E*cos(t*omega/(2*n))**2 * sin(t*omega) + E*cos(t*omega)*cos(t*omega/(2*n))*sin(t*omega/(2*n))/n

def fcx(x, y):
    return -q*x/(x**2 + y**2 + ac**2)**1.5

def fcy(x, y):
    return -q*y/(x**2 + y**2 + ac**2)**1.5

def equations(W, t):
    px, py, x, y = W
    return array([
        -El(t) + fcx(x, y),# == px'
        fcy(x, y),         # == py'
        px,                # == x'
        py                 # == y'
    ])

def momentum(v0, t0):
    # t0 = t0d*2*pi/(360*omega)
    tsp = linspace(t0, n*pi/omega, 1000)
    
    #[px_0, py_0, x_0, y_0]
    W0 = [0, v0, -ionpot/El(t0), 0]
    soln = integrate.odeint(equations, W0, tsp)[-1]
    
    if (soln[0]**2 + soln[1]**2)/2 + 1/(soln[2]**2 + soln[3]**2)**0.5:
        return soln
    else:
        return [0,0,0,0]

def nearest_id(array, value):
    return abs(array-value).argmin()

def tunneling_density(t0,v0_signed):
    f = abs(El(t0))
    F_a = (2*ionpot)**1.5
    W = 8*(2*ionpot)**0.5*F_a/f*exp(-2/3.0*F_a/f)
    a = (f/F_a)**0.5*(2*ionpot)**0.5
    w_perp = exp(-1*v0_signed**2/a**2)/(pi*a**2)
    return W*w_perp*pi

def generate_electron():
    while True:
        t0 = (t_end - t_start)*rand() + t_start
        v0 = (v0_max - v0_min)*rand() + v0_min
        xi = rand()
        supremum = 0.01 #00113 #0.000112430432247
        if xi <= tunneling_density(t0, v0)/supremum:
            return [t0, v0]

# t0, v0 = generate_electron()

# while True:
# for i in range(1000):
#     px_e, py_e, x_e, y_e = momentum(v0,t0)
#     print(t0, v0, px_e, py_e, x_e, y_e)

print(momentum(0.01,80))







