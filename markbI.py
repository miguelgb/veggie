#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from veggie import PLOT

X = np.linspace(0,10,256)
XY = PLOT(X, np.sin(X))
XY.Plot()

 
