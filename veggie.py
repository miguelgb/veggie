#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

class PLOT:
 def __init__(self, X, Y):
  self.X = X
  self.Y = Y

 def Plot(self):
  
  pl = plt.plot(self.X,self.Y)
  mng = plt.get_current_fig_manager() 
  mng.full_screen_toggle()
  plt.show()
  return pl

class metnum:
 """Numerical methods list
 4th Runge Kutt ----> RK4
 Example:
 fx1===>xi derivative function with arguments: (t, x, vx)
 fx2===>vi derivative function with arguments: (t, x, vx)
 f = [fx1, fx2]
 var = [time, xi, vxi]
 a = metnum(h, var, f)
 [time, xi, vxi] = a.RK42D()"""
 def __init__(self, h, x, f):
  self.h = h
  self.x = x
  self.f = f

 def RK41D(self):

  a1 = self.f(self.x[0],self.x[1])*self.h
  xk0 = self.x[0] + 0.5*self.h
  xk1 = self.x[1] + 0.5*a1

  a2 = self.f(xk0, xk1)*self.h
  xk0 = self.x[0] + 0.5*self.h
  xk1 = self.x[1] + 0.5*a2

  a3 = self.f(xk0, xk1)*self.h
  xk0 = self.x[0] + self.h
  xk1 = self.x[1] + a3

  a4 = self.f(xk0, xk1)*self.h
  self.x[0] = self.x[0] + self.h
  self.x[1] = self.x[1] + (a1 + 2*(a2 + a3) + a4)/6 
  return self.x 


 def RK42D(self):

  a1 = self.f[0](self.x[0],self.x[1],self.x[2])*self.h
  b1 = self.f[1](self.x[0],self.x[1],self.x[2])*self.h
  xk0 = self.x[0] + 0.5*self.h
  xk1 = self.x[1] + 0.5*a1
  xk2 = self.x[2] + 0.5*b1

  a2 = self.f[0](xk0, xk1, xk2)*self.h
  b2 = self.f[1](xk0, xk1, xk2)*self.h
  xk0 = self.x[0] + 0.5*self.h
  xk1 = self.x[1] + 0.5*a2
  xk2 = self.x[2] + 0.5*b2

  a3 = self.f[0](xk0, xk1, xk2)*self.h
  b3 = self.f[1](xk0, xk1, xk2)*self.h
  xk0 = self.x[0] + self.h
  xk1 = self.x[1] + a3
  xk2 = self.x[2] + b3

  a4 = self.f[0](xk0, xk1, xk2)*self.h
  b4 = self.f[1](xk0, xk1, xk2)*self.h
  self.x[0] = self.x[0] + self.h
  self.x[1] = self.x[1] + (a1 + 2*(a2 + a3) + a4)/6 
  self.x[2] = self.x[2] + (b1 + 2*(b2 + b3) + b4)/6 
  return self.x 


 def RK43D(self):

  a1 = self.f[0](self.x[0],self.x[1],self.x[2],self.x[3])*self.h
  b1 = self.f[1](self.x[0],self.x[1],self.x[2],self.x[3])*self.h
  c1 = self.f[2](self.x[0],self.x[1],self.x[2],self.x[3])*self.h
  xk0 = self.x[0] + 0.5*self.h
  xk1 = self.x[1] + 0.5*a1
  xk2 = self.x[2] + 0.5*b1
  xk3 = self.x[3] + 0.5*c1

  a2 = self.f[0](xk0, xk1, xk2, xk3)*self.h
  b2 = self.f[1](xk0, xk1, xk2, xk3)*self.h
  c2 = self.f[2](xk0, xk1, xk2, xk3)*self.h
  xk0 = self.x[0] + 0.5*self.h
  xk1 = self.x[1] + 0.5*a2
  xk2 = self.x[2] + 0.5*b2
  xk3 = self.x[3] + 0.5*c2

  a3 = self.f[0](xk0, xk1, xk2, xk3)*self.h
  b3 = self.f[1](xk0, xk1, xk2, xk3)*self.h
  c3 = self.f[2](xk0, xk1, xk2, xk3)*self.h
  xk0 = self.x[0] + self.h
  xk1 = self.x[1] + a3
  xk2 = self.x[2] + b3
  xk3 = self.x[3] + c3

  a4 = self.f[0](xk0, xk1, xk2, xk3)*self.h
  b4 = self.f[1](xk0, xk1, xk2, xk3)*self.h
  c4 = self.f[2](xk0, xk1, xk2, xk3)*self.h
  self.x[0] = self.x[0] + self.h
  self.x[1] = self.x[1] + (a1 + 2*(a2 + a3) + a4)/6 
  self.x[2] = self.x[2] + (b1 + 2*(b2 + b3) + b4)/6 
  self.x[3] = self.x[3] + (c1 + 2*(c2 + c3) + c4)/6 
  return self.x 


 def RK44D(self):

  a1 = self.f[0](self.x[0],self.x[1],self.x[2],self.x[3],self.x[4])*self.h
  b1 = self.f[1](self.x[0],self.x[1],self.x[2],self.x[3],self.x[4])*self.h
  c1 = self.f[2](self.x[0],self.x[1],self.x[2],self.x[3],self.x[4])*self.h
  d1 = self.f[3](self.x[0],self.x[1],self.x[2],self.x[3],self.x[4])*self.h
  xk0 = self.x[0] + 0.5*self.h
  xk1 = self.x[1] + 0.5*a1
  xk2 = self.x[2] + 0.5*b1
  xk3 = self.x[3] + 0.5*c1
  xk4 = self.x[4] + 0.5*d1

  a2 = self.f[0](xk0, xk1, xk2, xk3, xk4)*self.h
  b2 = self.f[1](xk0, xk1, xk2, xk3, xk4)*self.h
  c2 = self.f[2](xk0, xk1, xk2, xk3, xk4)*self.h
  d2 = self.f[3](xk0, xk1, xk2, xk3, xk4)*self.h
  xk0 = self.x[0] + 0.5*self.h
  xk1 = self.x[1] + 0.5*a2
  xk2 = self.x[2] + 0.5*b2
  xk3 = self.x[3] + 0.5*c2
  xk4 = self.x[4] + 0.5*d2

  a3 = self.f[0](xk0, xk1, xk2, xk3, xk4)*self.h
  b3 = self.f[1](xk0, xk1, xk2, xk3, xk4)*self.h
  c3 = self.f[2](xk0, xk1, xk2, xk3, xk4)*self.h
  d3 = self.f[3](xk0, xk1, xk2, xk3, xk4)*self.h
  xk0 = self.x[0] + self.h
  xk1 = self.x[1] + a3
  xk2 = self.x[2] + b3
  xk3 = self.x[3] + c3
  xk4 = self.x[4] + d3

  a4 = self.f[0](xk0, xk1, xk2, xk3, xk4)*self.h
  b4 = self.f[1](xk0, xk1, xk2, xk3, xk4)*self.h
  c4 = self.f[2](xk0, xk1, xk2, xk3, xk4)*self.h
  d4 = self.f[3](xk0, xk1, xk2, xk3, xk4)*self.h
  self.x[0] = self.x[0] + self.h
  self.x[1] = self.x[1] + (a1 + 2*(a2 + a3) + a4)/6 
  self.x[2] = self.x[2] + (b1 + 2*(b2 + b3) + b4)/6 
  self.x[3] = self.x[3] + (c1 + 2*(c2 + c3) + c4)/6 
  self.x[4] = self.x[4] + (d1 + 2*(d2 + d3) + d4)/6 
  return self.x 


 def RK45D(self):

  a1 = self.f[0](self.x[0],self.x[1],self.x[2],self.x[3],self.x[4],self.x[5])*self.h
  b1 = self.f[1](self.x[0],self.x[1],self.x[2],self.x[3],self.x[4],self.x[5])*self.h
  c1 = self.f[2](self.x[0],self.x[1],self.x[2],self.x[3],self.x[4],self.x[5])*self.h
  d1 = self.f[3](self.x[0],self.x[1],self.x[2],self.x[3],self.x[4],self.x[5])*self.h
  e1 = self.f[4](self.x[0],self.x[1],self.x[2],self.x[3],self.x[4],self.x[5])*self.h
  xk0 = self.x[0] + 0.5*self.h
  xk1 = self.x[1] + 0.5*a1
  xk2 = self.x[2] + 0.5*b1
  xk3 = self.x[3] + 0.5*c1
  xk4 = self.x[4] + 0.5*d1
  xk5 = self.x[5] + 0.5*e1

  a2 = self.f[0](xk0, xk1, xk2, xk3, xk4, xk5)*self.h
  b2 = self.f[1](xk0, xk1, xk2, xk3, xk4, xk5)*self.h
  c2 = self.f[2](xk0, xk1, xk2, xk3, xk4, xk5)*self.h
  d2 = self.f[3](xk0, xk1, xk2, xk3, xk4, xk5)*self.h
  e2 = self.f[4](xk0, xk1, xk2, xk3, xk4, xk5)*self.h
  xk0 = self.x[0] + 0.5*self.h
  xk1 = self.x[1] + 0.5*a2
  xk2 = self.x[2] + 0.5*b2
  xk3 = self.x[3] + 0.5*c2
  xk4 = self.x[4] + 0.5*d2
  xk5 = self.x[5] + 0.5*e2

  a3 = self.f[0](xk0, xk1, xk2, xk3, xk4, xk5)*self.h
  b3 = self.f[1](xk0, xk1, xk2, xk3, xk4, xk5)*self.h
  c3 = self.f[2](xk0, xk1, xk2, xk3, xk4, xk5)*self.h
  d3 = self.f[3](xk0, xk1, xk2, xk3, xk4, xk5)*self.h
  e3 = self.f[4](xk0, xk1, xk2, xk3, xk4, xk5)*self.h
  xk0 = self.x[0] + self.h
  xk1 = self.x[1] + a3
  xk2 = self.x[2] + b3
  xk3 = self.x[3] + c3
  xk4 = self.x[4] + d3
  xk5 = self.x[5] + e3

  a4 = self.f[0](xk0, xk1, xk2, xk3, xk4, xk5)*self.h
  b4 = self.f[1](xk0, xk1, xk2, xk3, xk4, xk5)*self.h
  c4 = self.f[2](xk0, xk1, xk2, xk3, xk4, xk5)*self.h
  d4 = self.f[3](xk0, xk1, xk2, xk3, xk4, xk5)*self.h
  e4 = self.f[4](xk0, xk1, xk2, xk3, xk4, xk5)*self.h
  self.x[0] = self.x[0] + self.h
  self.x[1] = self.x[1] + (a1 + 2*(a2 + a3) + a4)/6 
  self.x[2] = self.x[2] + (b1 + 2*(b2 + b3) + b4)/6 
  self.x[3] = self.x[3] + (c1 + 2*(c2 + c3) + c4)/6 
  self.x[4] = self.x[4] + (d1 + 2*(d2 + d3) + d4)/6 
  self.x[5] = self.x[5] + (e1 + 2*(e2 + e3) + e4)/6 
  return self.x 


 def RK46D(self):

  a1 = self.f[0](self.x[0],self.x[1],self.x[2],self.x[3],self.x[4],self.x[5],self.x[6])*self.h
  b1 = self.f[1](self.x[0],self.x[1],self.x[2],self.x[3],self.x[4],self.x[5],self.x[6])*self.h
  c1 = self.f[2](self.x[0],self.x[1],self.x[2],self.x[3],self.x[4],self.x[5],self.x[6])*self.h
  d1 = self.f[3](self.x[0],self.x[1],self.x[2],self.x[3],self.x[4],self.x[5],self.x[6])*self.h
  e1 = self.f[4](self.x[0],self.x[1],self.x[2],self.x[3],self.x[4],self.x[5],self.x[6])*self.h
  f1 = self.f[5](self.x[0],self.x[1],self.x[2],self.x[3],self.x[4],self.x[5],self.x[6])*self.h
  xk0 = self.x[0] + 0.5*self.h
  xk1 = self.x[1] + 0.5*a1
  xk2 = self.x[2] + 0.5*b1
  xk3 = self.x[3] + 0.5*c1
  xk4 = self.x[4] + 0.5*d1
  xk5 = self.x[5] + 0.5*e1
  xk6 = self.x[6] + 0.5*f1

  a2 = self.f[0](xk0, xk1, xk2, xk3, xk4, xk5, xk6)*self.h
  b2 = self.f[1](xk0, xk1, xk2, xk3, xk4, xk5, xk6)*self.h
  c2 = self.f[2](xk0, xk1, xk2, xk3, xk4, xk5, xk6)*self.h
  d2 = self.f[3](xk0, xk1, xk2, xk3, xk4, xk5, xk6)*self.h
  e2 = self.f[4](xk0, xk1, xk2, xk3, xk4, xk5, xk6)*self.h
  f2 = self.f[5](xk0, xk1, xk2, xk3, xk4, xk5, xk6)*self.h
  xk0 = self.x[0] + 0.5*self.h
  xk1 = self.x[1] + 0.5*a2
  xk2 = self.x[2] + 0.5*b2
  xk3 = self.x[3] + 0.5*c2
  xk4 = self.x[4] + 0.5*d2
  xk5 = self.x[5] + 0.5*e2
  xk6 = self.x[6] + 0.5*f2

  a3 = self.f[0](xk0, xk1, xk2, xk3, xk4, xk5, xk6)*self.h
  b3 = self.f[1](xk0, xk1, xk2, xk3, xk4, xk5, xk6)*self.h
  c3 = self.f[2](xk0, xk1, xk2, xk3, xk4, xk5, xk6)*self.h
  d3 = self.f[3](xk0, xk1, xk2, xk3, xk4, xk5, xk6)*self.h
  e3 = self.f[4](xk0, xk1, xk2, xk3, xk4, xk5, xk6)*self.h
  f3 = self.f[5](xk0, xk1, xk2, xk3, xk4, xk5, xk6)*self.h
  xk0 = self.x[0] + self.h
  xk1 = self.x[1] + a3
  xk2 = self.x[2] + b3
  xk3 = self.x[3] + c3
  xk4 = self.x[4] + d3
  xk5 = self.x[5] + e3
  xk6 = self.x[6] + f3

  a4 = self.f[0](xk0, xk1, xk2, xk3, xk4, xk5, xk6)*self.h
  b4 = self.f[1](xk0, xk1, xk2, xk3, xk4, xk5, xk6)*self.h
  c4 = self.f[2](xk0, xk1, xk2, xk3, xk4, xk5, xk6)*self.h
  d4 = self.f[3](xk0, xk1, xk2, xk3, xk4, xk5, xk6)*self.h
  e4 = self.f[4](xk0, xk1, xk2, xk3, xk4, xk5, xk6)*self.h
  f4 = self.f[5](xk0, xk1, xk2, xk3, xk4, xk5, xk6)*self.h
  self.x[0] = self.x[0] + self.h
  self.x[1] = self.x[1] + (a1 + 2*(a2 + a3) + a4)/6 
  self.x[2] = self.x[2] + (b1 + 2*(b2 + b3) + b4)/6 
  self.x[3] = self.x[3] + (c1 + 2*(c2 + c3) + c4)/6 
  self.x[4] = self.x[4] + (d1 + 2*(d2 + d3) + d4)/6 
  self.x[5] = self.x[5] + (e1 + 2*(e2 + e3) + e4)/6 
  self.x[6] = self.x[6] + (f1 + 2*(f2 + f3) + f4)/6 
  return self.x 


 def RK47D(self):

  a1 = self.f[0](self.x[0],self.x[1],self.x[2],self.x[3],self.x[4],self.x[5],self.x[6],self.x[7])*self.h
  b1 = self.f[1](self.x[0],self.x[1],self.x[2],self.x[3],self.x[4],self.x[5],self.x[6],self.x[7])*self.h
  c1 = self.f[2](self.x[0],self.x[1],self.x[2],self.x[3],self.x[4],self.x[5],self.x[6],self.x[7])*self.h
  d1 = self.f[3](self.x[0],self.x[1],self.x[2],self.x[3],self.x[4],self.x[5],self.x[6],self.x[7])*self.h
  e1 = self.f[4](self.x[0],self.x[1],self.x[2],self.x[3],self.x[4],self.x[5],self.x[6],self.x[7])*self.h
  f1 = self.f[5](self.x[0],self.x[1],self.x[2],self.x[3],self.x[4],self.x[5],self.x[6],self.x[7])*self.h
  g1 = self.f[6](self.x[0],self.x[1],self.x[2],self.x[3],self.x[4],self.x[5],self.x[6],self.x[7])*self.h
  xk0 = self.x[0] + 0.5*self.h
  xk1 = self.x[1] + 0.5*a1
  xk2 = self.x[2] + 0.5*b1
  xk3 = self.x[3] + 0.5*c1
  xk4 = self.x[4] + 0.5*d1
  xk5 = self.x[5] + 0.5*e1
  xk6 = self.x[6] + 0.5*f1
  xk7 = self.x[7] + 0.5*g1

  a2 = self.f[0](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7)*self.h
  b2 = self.f[1](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7)*self.h
  c2 = self.f[2](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7)*self.h
  d2 = self.f[3](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7)*self.h
  e2 = self.f[4](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7)*self.h
  f2 = self.f[5](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7)*self.h
  g2 = self.f[6](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7)*self.h
  xk0 = self.x[0] + 0.5*self.h
  xk1 = self.x[1] + 0.5*a2
  xk2 = self.x[2] + 0.5*b2
  xk3 = self.x[3] + 0.5*c2
  xk4 = self.x[4] + 0.5*d2
  xk5 = self.x[5] + 0.5*e2
  xk6 = self.x[6] + 0.5*f2
  xk7 = self.x[7] + 0.5*g2

  a3 = self.f[0](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7)*self.h
  b3 = self.f[1](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7)*self.h
  c3 = self.f[2](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7)*self.h
  d3 = self.f[3](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7)*self.h
  e3 = self.f[4](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7)*self.h
  f3 = self.f[5](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7)*self.h
  g3 = self.f[6](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7)*self.h
  xk0 = self.x[0] + self.h
  xk1 = self.x[1] + a3
  xk2 = self.x[2] + b3
  xk3 = self.x[3] + c3
  xk4 = self.x[4] + d3
  xk5 = self.x[5] + e3
  xk6 = self.x[6] + f3
  xk7 = self.x[7] + g3

  a4 = self.f[0](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7)*self.h
  b4 = self.f[1](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7)*self.h
  c4 = self.f[2](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7)*self.h
  d4 = self.f[3](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7)*self.h
  e4 = self.f[4](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7)*self.h
  f4 = self.f[5](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7)*self.h
  g4 = self.f[6](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7)*self.h
  self.x[0] = self.x[0] + self.h
  self.x[1] = self.x[1] + (a1 + 2*(a2 + a3) + a4)/6 
  self.x[2] = self.x[2] + (b1 + 2*(b2 + b3) + b4)/6 
  self.x[3] = self.x[3] + (c1 + 2*(c2 + c3) + c4)/6 
  self.x[4] = self.x[4] + (d1 + 2*(d2 + d3) + d4)/6 
  self.x[5] = self.x[5] + (e1 + 2*(e2 + e3) + e4)/6 
  self.x[6] = self.x[6] + (f1 + 2*(f2 + f3) + f4)/6 
  self.x[7] = self.x[7] + (g1 + 2*(g2 + g3) + g4)/6 
  return self.x 


 def RK48D(self):

  a1 = self.f[0](self.x[0],self.x[1],self.x[2],self.x[3],self.x[4],self.x[5],self.x[6],self.x[7],self.x[8])*self.h
  b1 = self.f[1](self.x[0],self.x[1],self.x[2],self.x[3],self.x[4],self.x[5],self.x[6],self.x[7],self.x[8])*self.h
  c1 = self.f[2](self.x[0],self.x[1],self.x[2],self.x[3],self.x[4],self.x[5],self.x[6],self.x[7],self.x[8])*self.h
  d1 = self.f[3](self.x[0],self.x[1],self.x[2],self.x[3],self.x[4],self.x[5],self.x[6],self.x[7],self.x[8])*self.h
  e1 = self.f[4](self.x[0],self.x[1],self.x[2],self.x[3],self.x[4],self.x[5],self.x[6],self.x[7],self.x[8])*self.h
  f1 = self.f[5](self.x[0],self.x[1],self.x[2],self.x[3],self.x[4],self.x[5],self.x[6],self.x[7],self.x[8])*self.h
  g1 = self.f[6](self.x[0],self.x[1],self.x[2],self.x[3],self.x[4],self.x[5],self.x[6],self.x[7],self.x[8])*self.h
  h1 = self.f[7](self.x[0],self.x[1],self.x[2],self.x[3],self.x[4],self.x[5],self.x[6],self.x[7],self.x[8])*self.h
  xk0 = self.x[0] + 0.5*self.h
  xk1 = self.x[1] + 0.5*a1
  xk2 = self.x[2] + 0.5*b1
  xk3 = self.x[3] + 0.5*c1
  xk4 = self.x[4] + 0.5*d1
  xk5 = self.x[5] + 0.5*e1
  xk6 = self.x[6] + 0.5*f1
  xk7 = self.x[7] + 0.5*g1
  xk8 = self.x[8] + 0.5*h1

  a2 = self.f[0](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8)*self.h
  b2 = self.f[1](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8)*self.h
  c2 = self.f[2](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8)*self.h
  d2 = self.f[3](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8)*self.h
  e2 = self.f[4](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8)*self.h
  f2 = self.f[5](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8)*self.h
  g2 = self.f[6](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8)*self.h
  h2 = self.f[7](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8)*self.h
  xk0 = self.x[0] + 0.5*self.h
  xk1 = self.x[1] + 0.5*a2
  xk2 = self.x[2] + 0.5*b2
  xk3 = self.x[3] + 0.5*c2
  xk4 = self.x[4] + 0.5*d2
  xk5 = self.x[5] + 0.5*e2
  xk6 = self.x[6] + 0.5*f2
  xk7 = self.x[7] + 0.5*g2
  xk8 = self.x[8] + 0.5*h2

  a3 = self.f[0](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8)*self.h
  b3 = self.f[1](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8)*self.h
  c3 = self.f[2](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8)*self.h
  d3 = self.f[3](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8)*self.h
  e3 = self.f[4](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8)*self.h
  f3 = self.f[5](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8)*self.h
  g3 = self.f[6](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8)*self.h
  h3 = self.f[7](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8)*self.h
  xk0 = self.x[0] + self.h
  xk1 = self.x[1] + a3
  xk2 = self.x[2] + b3
  xk3 = self.x[3] + c3
  xk4 = self.x[4] + d3
  xk5 = self.x[5] + e3
  xk6 = self.x[6] + f3
  xk7 = self.x[7] + g3
  xk8 = self.x[8] + h3

  a4 = self.f[0](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8)*self.h
  b4 = self.f[1](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8)*self.h
  c4 = self.f[2](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8)*self.h
  d4 = self.f[3](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8)*self.h
  e4 = self.f[4](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8)*self.h
  f4 = self.f[5](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8)*self.h
  g4 = self.f[6](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8)*self.h
  h4 = self.f[7](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8)*self.h
  self.x[0] = self.x[0] + self.h
  self.x[1] = self.x[1] + (a1 + 2*(a2 + a3) + a4)/6 
  self.x[2] = self.x[2] + (b1 + 2*(b2 + b3) + b4)/6 
  self.x[3] = self.x[3] + (c1 + 2*(c2 + c3) + c4)/6 
  self.x[4] = self.x[4] + (d1 + 2*(d2 + d3) + d4)/6 
  self.x[5] = self.x[5] + (e1 + 2*(e2 + e3) + e4)/6 
  self.x[6] = self.x[6] + (f1 + 2*(f2 + f3) + f4)/6 
  self.x[7] = self.x[7] + (g1 + 2*(g2 + g3) + g4)/6 
  self.x[8] = self.x[8] + (h1 + 2*(h2 + h3) + h4)/6 
  return self.x 


 def RK49D(self):

  a1 = self.f[0](self.x[0],self.x[1],self.x[2],self.x[3],self.x[4],self.x[5],self.x[6],self.x[7],self.x[8],self.x[9])*self.h
  b1 = self.f[1](self.x[0],self.x[1],self.x[2],self.x[3],self.x[4],self.x[5],self.x[6],self.x[7],self.x[8],self.x[9])*self.h
  c1 = self.f[2](self.x[0],self.x[1],self.x[2],self.x[3],self.x[4],self.x[5],self.x[6],self.x[7],self.x[8],self.x[9])*self.h
  d1 = self.f[3](self.x[0],self.x[1],self.x[2],self.x[3],self.x[4],self.x[5],self.x[6],self.x[7],self.x[8],self.x[9])*self.h
  e1 = self.f[4](self.x[0],self.x[1],self.x[2],self.x[3],self.x[4],self.x[5],self.x[6],self.x[7],self.x[8],self.x[9])*self.h
  f1 = self.f[5](self.x[0],self.x[1],self.x[2],self.x[3],self.x[4],self.x[5],self.x[6],self.x[7],self.x[8],self.x[9])*self.h
  g1 = self.f[6](self.x[0],self.x[1],self.x[2],self.x[3],self.x[4],self.x[5],self.x[6],self.x[7],self.x[8],self.x[9])*self.h
  h1 = self.f[7](self.x[0],self.x[1],self.x[2],self.x[3],self.x[4],self.x[5],self.x[6],self.x[7],self.x[8],self.x[9])*self.h
  i1 = self.f[8](self.x[0],self.x[1],self.x[2],self.x[3],self.x[4],self.x[5],self.x[6],self.x[7],self.x[8],self.x[9])*self.h
  xk0 = self.x[0] + 0.5*self.h
  xk1 = self.x[1] + 0.5*a1
  xk2 = self.x[2] + 0.5*b1
  xk3 = self.x[3] + 0.5*c1
  xk4 = self.x[4] + 0.5*d1
  xk5 = self.x[5] + 0.5*e1
  xk6 = self.x[6] + 0.5*f1
  xk7 = self.x[7] + 0.5*g1
  xk8 = self.x[8] + 0.5*h1
  xk9 = self.x[9] + 0.5*i1

  a2 = self.f[0](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9)*self.h
  b2 = self.f[1](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9)*self.h
  c2 = self.f[2](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9)*self.h
  d2 = self.f[3](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9)*self.h
  e2 = self.f[4](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9)*self.h
  f2 = self.f[5](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9)*self.h
  g2 = self.f[6](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9)*self.h
  h2 = self.f[7](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9)*self.h
  i2 = self.f[8](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9)*self.h
  xk0 = self.x[0] + 0.5*self.h
  xk1 = self.x[1] + 0.5*a2
  xk2 = self.x[2] + 0.5*b2
  xk3 = self.x[3] + 0.5*c2
  xk4 = self.x[4] + 0.5*d2
  xk5 = self.x[5] + 0.5*e2
  xk6 = self.x[6] + 0.5*f2
  xk7 = self.x[7] + 0.5*g2
  xk8 = self.x[8] + 0.5*h2
  xk9 = self.x[9] + 0.5*i2

  a3 = self.f[0](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9)*self.h
  b3 = self.f[1](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9)*self.h
  c3 = self.f[2](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9)*self.h
  d3 = self.f[3](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9)*self.h
  e3 = self.f[4](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9)*self.h
  f3 = self.f[5](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9)*self.h
  g3 = self.f[6](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9)*self.h
  h3 = self.f[7](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9)*self.h
  i3 = self.f[8](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9)*self.h
  xk0 = self.x[0] + self.h
  xk1 = self.x[1] + a3
  xk2 = self.x[2] + b3
  xk3 = self.x[3] + c3
  xk4 = self.x[4] + d3
  xk5 = self.x[5] + e3
  xk6 = self.x[6] + f3
  xk7 = self.x[7] + g3
  xk8 = self.x[8] + h3
  xk9 = self.x[9] + i3

  a4 = self.f[0](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9)*self.h
  b4 = self.f[1](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9)*self.h
  c4 = self.f[2](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9)*self.h
  d4 = self.f[3](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9)*self.h
  e4 = self.f[4](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9)*self.h
  f4 = self.f[5](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9)*self.h
  g4 = self.f[6](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9)*self.h
  h4 = self.f[7](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9)*self.h
  i4 = self.f[8](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9)*self.h
  self.x[0] = self.x[0] + self.h
  self.x[1] = self.x[1] + (a1 + 2*(a2 + a3) + a4)/6 
  self.x[2] = self.x[2] + (b1 + 2*(b2 + b3) + b4)/6 
  self.x[3] = self.x[3] + (c1 + 2*(c2 + c3) + c4)/6 
  self.x[4] = self.x[4] + (d1 + 2*(d2 + d3) + d4)/6 
  self.x[5] = self.x[5] + (e1 + 2*(e2 + e3) + e4)/6 
  self.x[6] = self.x[6] + (f1 + 2*(f2 + f3) + f4)/6 
  self.x[7] = self.x[7] + (g1 + 2*(g2 + g3) + g4)/6 
  self.x[8] = self.x[8] + (h1 + 2*(h2 + h3) + h4)/6 
  self.x[9] = self.x[9] + (i1 + 2*(i2 + i3) + i4)/6 
  return self.x 


 def RK410D(self):

  a1 = self.f[0](self.x[0],self.x[1],self.x[2],self.x[3],self.x[4],self.x[5],self.x[6],self.x[7],self.x[8],self.x[9],self.x[10])*self.h
  b1 = self.f[1](self.x[0],self.x[1],self.x[2],self.x[3],self.x[4],self.x[5],self.x[6],self.x[7],self.x[8],self.x[9],self.x[10])*self.h
  c1 = self.f[2](self.x[0],self.x[1],self.x[2],self.x[3],self.x[4],self.x[5],self.x[6],self.x[7],self.x[8],self.x[9],self.x[10])*self.h
  d1 = self.f[3](self.x[0],self.x[1],self.x[2],self.x[3],self.x[4],self.x[5],self.x[6],self.x[7],self.x[8],self.x[9],self.x[10])*self.h
  e1 = self.f[4](self.x[0],self.x[1],self.x[2],self.x[3],self.x[4],self.x[5],self.x[6],self.x[7],self.x[8],self.x[9],self.x[10])*self.h
  f1 = self.f[5](self.x[0],self.x[1],self.x[2],self.x[3],self.x[4],self.x[5],self.x[6],self.x[7],self.x[8],self.x[9],self.x[10])*self.h
  g1 = self.f[6](self.x[0],self.x[1],self.x[2],self.x[3],self.x[4],self.x[5],self.x[6],self.x[7],self.x[8],self.x[9],self.x[10])*self.h
  h1 = self.f[7](self.x[0],self.x[1],self.x[2],self.x[3],self.x[4],self.x[5],self.x[6],self.x[7],self.x[8],self.x[9],self.x[10])*self.h
  i1 = self.f[8](self.x[0],self.x[1],self.x[2],self.x[3],self.x[4],self.x[5],self.x[6],self.x[7],self.x[8],self.x[9],self.x[10])*self.h
  j1 = self.f[9](self.x[0],self.x[1],self.x[2],self.x[3],self.x[4],self.x[5],self.x[6],self.x[7],self.x[8],self.x[9],self.x[10])*self.h
  xk0 = self.x[0] + 0.5*self.h
  xk1 = self.x[1] + 0.5*a1
  xk2 = self.x[2] + 0.5*b1
  xk3 = self.x[3] + 0.5*c1
  xk4 = self.x[4] + 0.5*d1
  xk5 = self.x[5] + 0.5*e1
  xk6 = self.x[6] + 0.5*f1
  xk7 = self.x[7] + 0.5*g1
  xk8 = self.x[8] + 0.5*h1
  xk9 = self.x[9] + 0.5*i1
  xk10 = self.x[10] + 0.5*j1

  a2 = self.f[0](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10)*self.h
  b2 = self.f[1](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10)*self.h
  c2 = self.f[2](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10)*self.h
  d2 = self.f[3](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10)*self.h
  e2 = self.f[4](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10)*self.h
  f2 = self.f[5](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10)*self.h
  g2 = self.f[6](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10)*self.h
  h2 = self.f[7](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10)*self.h
  i2 = self.f[8](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10)*self.h
  j2 = self.f[9](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10)*self.h
  xk0 = self.x[0] + 0.5*self.h
  xk1 = self.x[1] + 0.5*a2
  xk2 = self.x[2] + 0.5*b2
  xk3 = self.x[3] + 0.5*c2
  xk4 = self.x[4] + 0.5*d2
  xk5 = self.x[5] + 0.5*e2
  xk6 = self.x[6] + 0.5*f2
  xk7 = self.x[7] + 0.5*g2
  xk8 = self.x[8] + 0.5*h2
  xk9 = self.x[9] + 0.5*i2
  xk10 = self.x[10] + 0.5*j2

  a3 = self.f[0](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10)*self.h
  b3 = self.f[1](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10)*self.h
  c3 = self.f[2](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10)*self.h
  d3 = self.f[3](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10)*self.h
  e3 = self.f[4](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10)*self.h
  f3 = self.f[5](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10)*self.h
  g3 = self.f[6](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10)*self.h
  h3 = self.f[7](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10)*self.h
  i3 = self.f[8](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10)*self.h
  j3 = self.f[9](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10)*self.h
  xk0 = self.x[0] + self.h
  xk1 = self.x[1] + a3
  xk2 = self.x[2] + b3
  xk3 = self.x[3] + c3
  xk4 = self.x[4] + d3
  xk5 = self.x[5] + e3
  xk6 = self.x[6] + f3
  xk7 = self.x[7] + g3
  xk8 = self.x[8] + h3
  xk9 = self.x[9] + i3
  xk10 = self.x[10] + j3

  a4 = self.f[0](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10)*self.h
  b4 = self.f[1](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10)*self.h
  c4 = self.f[2](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10)*self.h
  d4 = self.f[3](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10)*self.h
  e4 = self.f[4](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10)*self.h
  f4 = self.f[5](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10)*self.h
  g4 = self.f[6](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10)*self.h
  h4 = self.f[7](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10)*self.h
  i4 = self.f[8](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10)*self.h
  j4 = self.f[9](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10)*self.h
  self.x[0] = self.x[0] + self.h
  self.x[1] = self.x[1] + (a1 + 2*(a2 + a3) + a4)/6 
  self.x[2] = self.x[2] + (b1 + 2*(b2 + b3) + b4)/6 
  self.x[3] = self.x[3] + (c1 + 2*(c2 + c3) + c4)/6 
  self.x[4] = self.x[4] + (d1 + 2*(d2 + d3) + d4)/6 
  self.x[5] = self.x[5] + (e1 + 2*(e2 + e3) + e4)/6 
  self.x[6] = self.x[6] + (f1 + 2*(f2 + f3) + f4)/6 
  self.x[7] = self.x[7] + (g1 + 2*(g2 + g3) + g4)/6 
  self.x[8] = self.x[8] + (h1 + 2*(h2 + h3) + h4)/6 
  self.x[9] = self.x[9] + (i1 + 2*(i2 + i3) + i4)/6 
  self.x[10] = self.x[10] + (j1 + 2*(j2 + j3) + j4)/6 
  return self.x 


 def RK411D(self):

  a1 = self.f[0](self.x[0],self.x[1],self.x[2],self.x[3],self.x[4],self.x[5],self.x[6],self.x[7],self.x[8],self.x[9],self.x[10],self.x[11])*self.h
  b1 = self.f[1](self.x[0],self.x[1],self.x[2],self.x[3],self.x[4],self.x[5],self.x[6],self.x[7],self.x[8],self.x[9],self.x[10],self.x[11])*self.h
  c1 = self.f[2](self.x[0],self.x[1],self.x[2],self.x[3],self.x[4],self.x[5],self.x[6],self.x[7],self.x[8],self.x[9],self.x[10],self.x[11])*self.h
  d1 = self.f[3](self.x[0],self.x[1],self.x[2],self.x[3],self.x[4],self.x[5],self.x[6],self.x[7],self.x[8],self.x[9],self.x[10],self.x[11])*self.h
  e1 = self.f[4](self.x[0],self.x[1],self.x[2],self.x[3],self.x[4],self.x[5],self.x[6],self.x[7],self.x[8],self.x[9],self.x[10],self.x[11])*self.h
  f1 = self.f[5](self.x[0],self.x[1],self.x[2],self.x[3],self.x[4],self.x[5],self.x[6],self.x[7],self.x[8],self.x[9],self.x[10],self.x[11])*self.h
  g1 = self.f[6](self.x[0],self.x[1],self.x[2],self.x[3],self.x[4],self.x[5],self.x[6],self.x[7],self.x[8],self.x[9],self.x[10],self.x[11])*self.h
  h1 = self.f[7](self.x[0],self.x[1],self.x[2],self.x[3],self.x[4],self.x[5],self.x[6],self.x[7],self.x[8],self.x[9],self.x[10],self.x[11])*self.h
  i1 = self.f[8](self.x[0],self.x[1],self.x[2],self.x[3],self.x[4],self.x[5],self.x[6],self.x[7],self.x[8],self.x[9],self.x[10],self.x[11])*self.h
  j1 = self.f[9](self.x[0],self.x[1],self.x[2],self.x[3],self.x[4],self.x[5],self.x[6],self.x[7],self.x[8],self.x[9],self.x[10],self.x[11])*self.h
  k1 = self.f[10](self.x[0],self.x[1],self.x[2],self.x[3],self.x[4],self.x[5],self.x[6],self.x[7],self.x[8],self.x[9],self.x[10],self.x[11])*self.h
  xk0 = self.x[0] + 0.5*self.h
  xk1 = self.x[1] + 0.5*a1
  xk2 = self.x[2] + 0.5*b1
  xk3 = self.x[3] + 0.5*c1
  xk4 = self.x[4] + 0.5*d1
  xk5 = self.x[5] + 0.5*e1
  xk6 = self.x[6] + 0.5*f1
  xk7 = self.x[7] + 0.5*g1
  xk8 = self.x[8] + 0.5*h1
  xk9 = self.x[9] + 0.5*i1
  xk10 = self.x[10] + 0.5*j1
  xk11 = self.x[11] + 0.5*k1

  a2 = self.f[0](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10, xk11)*self.h
  b2 = self.f[1](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10, xk11)*self.h
  c2 = self.f[2](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10, xk11)*self.h
  d2 = self.f[3](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10, xk11)*self.h
  e2 = self.f[4](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10, xk11)*self.h
  f2 = self.f[5](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10, xk11)*self.h
  g2 = self.f[6](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10, xk11)*self.h
  h2 = self.f[7](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10, xk11)*self.h
  i2 = self.f[8](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10, xk11)*self.h
  j2 = self.f[9](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10, xk11)*self.h
  k2 = self.f[10](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10, xk11)*self.h
  xk0 = self.x[0] + 0.5*self.h
  xk1 = self.x[1] + 0.5*a2
  xk2 = self.x[2] + 0.5*b2
  xk3 = self.x[3] + 0.5*c2
  xk4 = self.x[4] + 0.5*d2
  xk5 = self.x[5] + 0.5*e2
  xk6 = self.x[6] + 0.5*f2
  xk7 = self.x[7] + 0.5*g2
  xk8 = self.x[8] + 0.5*h2
  xk9 = self.x[9] + 0.5*i2
  xk10 = self.x[10] + 0.5*j2
  xk11 = self.x[11] + 0.5*k2

  a3 = self.f[0](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10, xk11)*self.h
  b3 = self.f[1](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10, xk11)*self.h
  c3 = self.f[2](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10, xk11)*self.h
  d3 = self.f[3](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10, xk11)*self.h
  e3 = self.f[4](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10, xk11)*self.h
  f3 = self.f[5](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10, xk11)*self.h
  g3 = self.f[6](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10, xk11)*self.h
  h3 = self.f[7](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10, xk11)*self.h
  i3 = self.f[8](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10, xk11)*self.h
  j3 = self.f[9](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10, xk11)*self.h
  k3 = self.f[10](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10, xk11)*self.h
  xk0 = self.x[0] + self.h
  xk1 = self.x[1] + a3
  xk2 = self.x[2] + b3
  xk3 = self.x[3] + c3
  xk4 = self.x[4] + d3
  xk5 = self.x[5] + e3
  xk6 = self.x[6] + f3
  xk7 = self.x[7] + g3
  xk8 = self.x[8] + h3
  xk9 = self.x[9] + i3
  xk10 = self.x[10] + j3
  xk11 = self.x[11] + k3

  a4 = self.f[0](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10, xk11)*self.h
  b4 = self.f[1](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10, xk11)*self.h
  c4 = self.f[2](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10, xk11)*self.h
  d4 = self.f[3](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10, xk11)*self.h
  e4 = self.f[4](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10, xk11)*self.h
  f4 = self.f[5](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10, xk11)*self.h
  g4 = self.f[6](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10, xk11)*self.h
  h4 = self.f[7](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10, xk11)*self.h
  i4 = self.f[8](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10, xk11)*self.h
  j4 = self.f[9](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10, xk11)*self.h
  k4 = self.f[10](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10, xk11)*self.h
  self.x[0] = self.x[0] + self.h
  self.x[1] = self.x[1] + (a1 + 2*(a2 + a3) + a4)/6 
  self.x[2] = self.x[2] + (b1 + 2*(b2 + b3) + b4)/6 
  self.x[3] = self.x[3] + (c1 + 2*(c2 + c3) + c4)/6 
  self.x[4] = self.x[4] + (d1 + 2*(d2 + d3) + d4)/6 
  self.x[5] = self.x[5] + (e1 + 2*(e2 + e3) + e4)/6 
  self.x[6] = self.x[6] + (f1 + 2*(f2 + f3) + f4)/6 
  self.x[7] = self.x[7] + (g1 + 2*(g2 + g3) + g4)/6 
  self.x[8] = self.x[8] + (h1 + 2*(h2 + h3) + h4)/6 
  self.x[9] = self.x[9] + (i1 + 2*(i2 + i3) + i4)/6 
  self.x[10] = self.x[10] + (j1 + 2*(j2 + j3) + j4)/6 
  self.x[11] = self.x[11] + (k1 + 2*(k2 + k3) + k4)/6 
  return self.x 


 def RK412D(self):

  a1 = self.f[0](self.x[0],self.x[1],self.x[2],self.x[3],self.x[4],self.x[5],self.x[6],self.x[7],self.x[8],self.x[9],self.x[10],self.x[11],self.x[12])*self.h
  b1 = self.f[1](self.x[0],self.x[1],self.x[2],self.x[3],self.x[4],self.x[5],self.x[6],self.x[7],self.x[8],self.x[9],self.x[10],self.x[11],self.x[12])*self.h
  c1 = self.f[2](self.x[0],self.x[1],self.x[2],self.x[3],self.x[4],self.x[5],self.x[6],self.x[7],self.x[8],self.x[9],self.x[10],self.x[11],self.x[12])*self.h
  d1 = self.f[3](self.x[0],self.x[1],self.x[2],self.x[3],self.x[4],self.x[5],self.x[6],self.x[7],self.x[8],self.x[9],self.x[10],self.x[11],self.x[12])*self.h
  e1 = self.f[4](self.x[0],self.x[1],self.x[2],self.x[3],self.x[4],self.x[5],self.x[6],self.x[7],self.x[8],self.x[9],self.x[10],self.x[11],self.x[12])*self.h
  f1 = self.f[5](self.x[0],self.x[1],self.x[2],self.x[3],self.x[4],self.x[5],self.x[6],self.x[7],self.x[8],self.x[9],self.x[10],self.x[11],self.x[12])*self.h
  g1 = self.f[6](self.x[0],self.x[1],self.x[2],self.x[3],self.x[4],self.x[5],self.x[6],self.x[7],self.x[8],self.x[9],self.x[10],self.x[11],self.x[12])*self.h
  h1 = self.f[7](self.x[0],self.x[1],self.x[2],self.x[3],self.x[4],self.x[5],self.x[6],self.x[7],self.x[8],self.x[9],self.x[10],self.x[11],self.x[12])*self.h
  i1 = self.f[8](self.x[0],self.x[1],self.x[2],self.x[3],self.x[4],self.x[5],self.x[6],self.x[7],self.x[8],self.x[9],self.x[10],self.x[11],self.x[12])*self.h
  j1 = self.f[9](self.x[0],self.x[1],self.x[2],self.x[3],self.x[4],self.x[5],self.x[6],self.x[7],self.x[8],self.x[9],self.x[10],self.x[11],self.x[12])*self.h
  k1 = self.f[10](self.x[0],self.x[1],self.x[2],self.x[3],self.x[4],self.x[5],self.x[6],self.x[7],self.x[8],self.x[9],self.x[10],self.x[11],self.x[12])*self.h
  l1 = self.f[11](self.x[0],self.x[1],self.x[2],self.x[3],self.x[4],self.x[5],self.x[6],self.x[7],self.x[8],self.x[9],self.x[10],self.x[11],self.x[12])*self.h
  xk0 = self.x[0] + 0.5*self.h
  xk1 = self.x[1] + 0.5*a1
  xk2 = self.x[2] + 0.5*b1
  xk3 = self.x[3] + 0.5*c1
  xk4 = self.x[4] + 0.5*d1
  xk5 = self.x[5] + 0.5*e1
  xk6 = self.x[6] + 0.5*f1
  xk7 = self.x[7] + 0.5*g1
  xk8 = self.x[8] + 0.5*h1
  xk9 = self.x[9] + 0.5*i1
  xk10 = self.x[10] + 0.5*j1
  xk11 = self.x[11] + 0.5*k1
  xk12 = self.x[12] + 0.5*l1

  a2 = self.f[0](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10, xk11, xk12)*self.h
  b2 = self.f[1](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10, xk11, xk12)*self.h
  c2 = self.f[2](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10, xk11, xk12)*self.h
  d2 = self.f[3](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10, xk11, xk12)*self.h
  e2 = self.f[4](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10, xk11, xk12)*self.h
  f2 = self.f[5](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10, xk11, xk12)*self.h
  g2 = self.f[6](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10, xk11, xk12)*self.h
  h2 = self.f[7](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10, xk11, xk12)*self.h
  i2 = self.f[8](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10, xk11, xk12)*self.h
  j2 = self.f[9](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10, xk11, xk12)*self.h
  k2 = self.f[10](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10, xk11, xk12)*self.h
  l2 = self.f[11](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10, xk11, xk12)*self.h
  xk0 = self.x[0] + 0.5*self.h
  xk1 = self.x[1] + 0.5*a2
  xk2 = self.x[2] + 0.5*b2
  xk3 = self.x[3] + 0.5*c2
  xk4 = self.x[4] + 0.5*d2
  xk5 = self.x[5] + 0.5*e2
  xk6 = self.x[6] + 0.5*f2
  xk7 = self.x[7] + 0.5*g2
  xk8 = self.x[8] + 0.5*h2
  xk9 = self.x[9] + 0.5*i2
  xk10 = self.x[10] + 0.5*j2
  xk11 = self.x[11] + 0.5*k2
  xk12 = self.x[12] + 0.5*l2

  a3 = self.f[0](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10, xk11, xk12)*self.h
  b3 = self.f[1](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10, xk11, xk12)*self.h
  c3 = self.f[2](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10, xk11, xk12)*self.h
  d3 = self.f[3](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10, xk11, xk12)*self.h
  e3 = self.f[4](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10, xk11, xk12)*self.h
  f3 = self.f[5](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10, xk11, xk12)*self.h
  g3 = self.f[6](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10, xk11, xk12)*self.h
  h3 = self.f[7](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10, xk11, xk12)*self.h
  i3 = self.f[8](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10, xk11, xk12)*self.h
  j3 = self.f[9](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10, xk11, xk12)*self.h
  k3 = self.f[10](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10, xk11, xk12)*self.h
  l3 = self.f[11](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10, xk11, xk12)*self.h
  xk0 = self.x[0] + self.h
  xk1 = self.x[1] + a3
  xk2 = self.x[2] + b3
  xk3 = self.x[3] + c3
  xk4 = self.x[4] + d3
  xk5 = self.x[5] + e3
  xk6 = self.x[6] + f3
  xk7 = self.x[7] + g3
  xk8 = self.x[8] + h3
  xk9 = self.x[9] + i3
  xk10 = self.x[10] + j3
  xk11 = self.x[11] + k3
  xk12 = self.x[12] + l3

  a4 = self.f[0](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10, xk11, xk12)*self.h
  b4 = self.f[1](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10, xk11, xk12)*self.h
  c4 = self.f[2](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10, xk11, xk12)*self.h
  d4 = self.f[3](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10, xk11, xk12)*self.h
  e4 = self.f[4](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10, xk11, xk12)*self.h
  f4 = self.f[5](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10, xk11, xk12)*self.h
  g4 = self.f[6](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10, xk11, xk12)*self.h
  h4 = self.f[7](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10, xk11, xk12)*self.h
  i4 = self.f[8](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10, xk11, xk12)*self.h
  j4 = self.f[9](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10, xk11, xk12)*self.h
  k4 = self.f[10](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10, xk11, xk12)*self.h
  l4 = self.f[11](xk0, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8, xk9, xk10, xk11, xk12)*self.h
  self.x[0] = self.x[0] + self.h
  self.x[1] = self.x[1] + (a1 + 2*(a2 + a3) + a4)/6 
  self.x[2] = self.x[2] + (b1 + 2*(b2 + b3) + b4)/6 
  self.x[3] = self.x[3] + (c1 + 2*(c2 + c3) + c4)/6 
  self.x[4] = self.x[4] + (d1 + 2*(d2 + d3) + d4)/6 
  self.x[5] = self.x[5] + (e1 + 2*(e2 + e3) + e4)/6 
  self.x[6] = self.x[6] + (f1 + 2*(f2 + f3) + f4)/6 
  self.x[7] = self.x[7] + (g1 + 2*(g2 + g3) + g4)/6 
  self.x[8] = self.x[8] + (h1 + 2*(h2 + h3) + h4)/6 
  self.x[9] = self.x[9] + (i1 + 2*(i2 + i3) + i4)/6 
  self.x[10] = self.x[10] + (j1 + 2*(j2 + j3) + j4)/6 
  self.x[11] = self.x[11] + (k1 + 2*(k2 + k3) + k4)/6 
  self.x[12] = self.x[12] + (l1 + 2*(l2 + l3) + l4)/6 
  return self.x 

