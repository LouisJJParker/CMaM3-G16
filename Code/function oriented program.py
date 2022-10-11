import math
from tkinter import Frame
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import start_conditions

L = start_conditions.L # string length
n_b = start_conditions.n_b # Number of oscillators
rho = start_conditions.rho # mass per unit length
alpha = start_conditions.alpha #non-linear coefficient
k_young = start_conditions.k_young # Young's modulus
init_type = start_conditions.init_type # initial conditions

dt = start_conditions.dt # time step
final_time = start_conditions.final_time # total simulation time

i = np.linspace(0,L,n_b)
x = np.sin(i*2*np.pi)
a = np.zeros_like(i)
v = np.zeros_like(i)

fig, ax = plt.subplots()

line, = ax.plot(i,x)

def simulate(a,v,x):
    for j in range(1,n_b-1):
        a[j] = ((k_young/rho) * 
                (x[j-1]+x[j+1]-2*x[j]) * 
                (1+alpha*(x[j+1]-x[j-1])))
    v += a*dt
    x += v*dt
    return x

def animate(t):
    y = simulate(a,v,x)
    line.set_ydata(y)
    return line,

def show_graph():
    ani = animation.FuncAnimation(
        fig, animate, interval=0.001, blit=True, save_count=50)

    plt.show()

show_graph()
