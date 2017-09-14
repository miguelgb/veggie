#!/usr/bin/env python 
# -*- coding: utf-8 -*-

# This code tests the 4th order Runge-Kutta integrator on veggie using a logaritmic
# bar potential

import math as mt
import numpy as np
import matplotlib.pyplot as plt
import time
import sys
from veggie import metnum #This calls the metnum library on veggie

# These are the initial values for differents trayectories 
vo2 = [1.0, 1.0, 1.0, 1.0, 3.0, 40.0, 1.0]
rc = [0.14, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0]
q = [0.9, 0.7, 0.9, 0.7, 1.0, 1.0, 0.7]
x0 = [1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0]
y0 = [1.3, 0.13, 0.13, 0.064, 0.0, 0.0, 0.2663]
vx0 = [0.0, 1.83496, 1.95121, 2.1575, 0.0, 0.0, 1.39029]
vy0 = [0.0, 0.0, 0.25, 0.36, 0.3, 6.28318, 0.0]

# When the script is run, an integer n-value between 0 and 6 must be written so the
# line on the terminal is "python markIII.py n" 
n = int(sys.argv[1])
h, tiny = 2.5e-5, 1e-18
e0 = 0.5*(vx0[n]**2 + vy0[n]**2) + 0.5*vo2[n]*mt.log(x0[n]**2 + (y0[n]/q[n])**2 + rc[n]**2)
T = int(100/h)
t = np.linspace(0, T, num=T+1)
xi, vxi, yi, vyi, ei, time = x0[n], vx0[n], y0[n], vy0[n], e0, 0.
X, Y, VX, VY, TIME, E = [], [], [], [], [], []
ages = 0.

def fx1(t, x, vx):
 return vx

def fy1(t, y, vy):
 return vy

def fx2(t, x, vx):
 return -(vo2[n]*x)/(x**2 + (yi/q[n])**2 + rc[n]**2)

def fy2(t, y, vy):
 return -(vo2[n]*y/q[n]**2)/(xi**2 + (y/q[n])**2 + rc[n]**2)

fx = [fx1, fx2]
fy = [fy1, fy2]

def error(x, vx, y, vy):
 err = 0.5*(vx**2 + vy**2) + 0.5*vo2[n]*mt.log(x**2 + (y/q[n])**2 + rc[n]**2)
 err = abs((err-e0)/e0)
 if err < tiny:
  err = -15
 else:
  err = mt.log(err,10)
 return err

for i in t:
 VX.append(vxi)
 X.append(xi)
 VY.append(vyi)
 Y.append(yi)
 TIME.append(time)
# varx/vary have [time, x/y-position value, x/y-velocity value]
 varx = [time, xi, vxi]
 vary = [time, yi, vyi] 
 ax = metnum(h, varx, fx)
 ay = metnum(h, vary, fy)
 [tim, xi, vxi] = ax.RK42D()
 [time, yi, vyi] = ay.RK42D()
 ei = error(xi, vxi, yi, vyi)
 E.append(ei)
 if ages<=time:
  print time,ei,e0
  ages=ages+1

# Plot the trayectory for a i given initial conditions
plt.figure(1)
plt.plot(X,Y)
plt.figure(2)
plt.plot(TIME,X)
# Plot the relative logarithmic error of the simulation
plt.figure(3)
plt.plot(TIME,E)
plt.grid()
plt.show()

