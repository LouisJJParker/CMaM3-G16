import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import start_conditions
from genericpath import exists

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

def main():
    #write_file(run_without_graph())
    graph_from_file("FPUT_experiment0.165.npy")
    #plot_at_t(10000, read_file("FPUT_experiment0.165.npy"))

def read_file(file_name):
    data = np.load(file = file_name)
    return data

def write_file(data):
    file_name = "FPUT_experiment{}{}.npy".format(dt,n_b)
    if exists(file_name):
        file=open(file_name,'w')
        file.close()
        print('x')
    np.save(file=file_name, arr=data)

def simulate(a,v,x):
    for j in range(1,n_b-1):
        a[j] = ((k_young/rho) * 
                (x[j-1]+x[j+1]-2*x[j]) * 
                (1+alpha*(x[j+1]-x[j-1])))
    v += a*dt
    x += v*dt
    return x

def run_without_graph():
    data = []
    for i in range(final_time):
        data.append(simulate(a,v,x).copy())
    return data

def animate(t):
    y = simulate(a,v,x)
    line.set_ydata(y)
    return line,

def animate_from_data(frame, data):
    y = data[frame]
    line.set_ydata(y)
    return line,

def show_graph():
    ani = animation.FuncAnimation(
        fig, animate, interval=0.001, blit=True, save_count=50)

    plt.show()

def graph_from_data(data):
    ani = animation.FuncAnimation(
        fig, animate_from_data, frames=range(0,final_time,40), fargs=[data], interval=0.1, blit=True, repeat = False)

    plt.show()    

def graph_from_file(file_name):
    data = read_file(file_name)

    ani = animation.FuncAnimation(
        fig, animate_from_data, frames=range(0,final_time), fargs=[data], interval=0.1, blit=True, repeat = False)

    plt.show()

def plot_at_t(t,data):
    line.set_ydata(data[t])

    plt.show()

main()