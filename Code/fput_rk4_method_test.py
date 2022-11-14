import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.animation import PillowWriter
import math

L = 1
N = 65
rho = 400
y_m = 0.1
alp = 0.25
amp = 1.0

T = 12000
dt = 0.1

n_step = math.ceil(T/dt)

h = L/(N-1)

m = rho*h

k = y_m/h

i = np.linspace(0, L, N)
x_i = np.zeros(len(i))
y2 = np.zeros(len(i))
y3 = np.zeros(len(i))
u_i = np.zeros(len(i))
a_i = np.zeros(len(i))

for n in range(N):
    x_i[n] = (amp)*np.sin(2*np.pi*i[n])
    
x_i[0] = x_i[-1] = 0

for n in range(1,N):
    y2[n] = x_i[n-1]

for n in range(1,N-1):
    y3[n] = x_i[n+1]

fig, ax = plt.subplots()
ln, = ax.plot(i, x_i)

#y1 = x_i[n]
#y2 = x_i[n-1]
#y3 = x_i[n+1]

def model(y1, y2, y3):
    dudt = (k/m) * (y2 + y3 -2*y1) * (1 + alp*(y3 - y2))
    return dudt

data = []
for t in range(n_step):
    
    for n in range(1,N-1):
        
        y_dum_1 = x_i[n]
        y_dum_2 = y2[n]
        y_dum_3 = y3[n]
        k1 = model(y_dum_1, y_dum_2, y_dum_3)
        
        y_dum_1 = x_i[n] + k1 * dt/2
        y_dum_2 = y2[n] + k1 * dt/2
        y_dum_3 = y3[n] + k1 * dt/2
        k2 = model(y_dum_1, y_dum_2, y_dum_3)
        
        y_dum_1 = x_i[n] + k2 * dt/2
        y_dum_2 = y2[n] + k2 * dt/2
        y_dum_3 = y3[n] + k2 * dt/2
        k3 = model(y_dum_1, y_dum_2, y_dum_3)
        
        y_dum_1 = x_i[n] + k3*dt
        y_dum_2 = y2[n] + k3*dt
        y_dum_3 = y3[n] + k3*dt
        k4 = model(y_dum_1, y_dum_2, y_dum_3)
        
        a_i[n] = 1/6 * k1 + 2/6 * k2 + 2/6 * k3 + 1/6 * k4
        
    u_i += a_i*dt
    x_i += u_i*dt
    for n in range(1,N):
        y2[n] = x_i[n-1]

    for n in range(1,N-1):
        y3[n] = x_i[n+1]
    
    data.append(x_i.copy())
    
def settings_data(fr, data):
    y = data[fr]
    ln.set_ydata(y)
    return ln,

def anime_graph_file(data):
    
    ani = FuncAnimation(
        fig, settings_data, frames=range(0,n_step,50), fargs=[data], interval=0.1, blit=True, repeat=False)
    
    sav = PillowWriter(fps=80)
    ani.save('Simulation_RK4_{}N.gif'.format(N), writer=sav)
    
    plt.show()

    
file_name = 'FPUT_RK4_Simulation_{}N.npy'.format(N)

np.save(file_name, arr=data)

anime_graph_file(data)
