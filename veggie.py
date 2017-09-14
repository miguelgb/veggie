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
