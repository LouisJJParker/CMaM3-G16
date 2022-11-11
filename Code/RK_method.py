import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import PillowWriter
from scipy.integrate import solve_ivp
import start_conditions
import os

L = start_conditions.L # string length
n_b = start_conditions.n_b+2 # Number of oscillators + the two fixed points
rho = start_conditions.rho # mass per unit length
alpha = start_conditions.alpha #non-linear coefficient
k_young = start_conditions.k_young # Young's modulus
init_type = start_conditions.init_type # initial conditions

dt = start_conditions.dt # time step
final_time = start_conditions.final_time # total simulation time in seconds
RK_n = start_conditions.RK_n # Number of slopes for RK solution

frame_count = int(final_time//dt) #number of simulated frames

h = L/(n_b-3) #maybe should be distance between x[i-1] and x[i+1]
i = np.linspace(0,L,n_b)
x = np.sin(i*2*np.pi)
v = np.zeros_like(i)
a = np.zeros_like(i)

data = []

fig, ax = plt.subplots()

line, = ax.plot(i,x)

def model(x): # Using Euler's method
    for j in range(1,n_b-1):
        a[j] = ((k_young/(rho*h**2)) * 
                (x[j-1]+x[j+1]-2*x[j]) * 
                (1+alpha*(x[j+1]-x[j-1])))
    return a

def graph_from_data(data):
    ani = animation.FuncAnimation(
        fig, animate_from_data, frames=int(frame_count), fargs=[data], interval=0.01, blit=True, repeat = False)

    writer = PillowWriter(fps=25)
    ani.save("demo_sine{}{}_RK.gif".format(dt,n_b), writer=writer)

    plt.show()

def animate_from_data(frame, data):
    y = data[frame]
    line.set_ydata(y)
    return line,

def write_file(data):
    file_name = "FPUT_experiment{}{}_RK.npy".format(dt,n_b)
    if os.path.exists(file_name):
        os.remove(file_name)
    np.save(file=file_name, arr=data)

for i in range(frame_count):
    data.append(x)

    y_dummy = x
    k1 = model(y_dummy)

    y_dummy = x + k1 * dt/2
    k2 = model(y_dummy)

    y_dummy = x + k2 * dt/2
    k3 = model(y_dummy)

    y_dummy = x + k3 * dt
    k4 = model(y_dummy)

    slope = 1/6*k1 + 2/6*k2 + 2/6*k3 + 1/6*k4

    v = v + dt*slope
    x = x + dt*v

write_file(data)
print('file saved')
graph_from_data(data)