from numpy import sin, cos
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

t_stop = 15
dt = 0.02
t = np.arange(0,t_stop,dt)

n_b =65

f = lambda t,x; 