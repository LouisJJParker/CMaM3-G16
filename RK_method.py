# Import essential built-in modules
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import PillowWriter
from scipy.integrate import solve_ivp
import start_conditions
import os


# Import initial conditions from file "start_conditions.py"
L = start_conditions.L # String length
n_b = start_conditions.n_b + 2 # Number of oscillators + the two fixed points
rho = start_conditions.rho # Mass per unit length
alpha = start_conditions.alpha # Non-linear coefficient
k_young = start_conditions.k_young # Young's modulus
init_type = start_conditions.init_type # Initial conditions

dt = start_conditions.dt # Time step
final_time = start_conditions.final_time # Total simulation time in seconds
RK_n = start_conditions.RK_n # Number of slopes for RK solution

frame_count = int(final_time//dt) # Number of simulated frames

h = L/(n_b - 3) # Latice spacing distance between two neighboring particles
i = np.linspace(0, L, n_b) # Indexes of particles
x = np.sin(i*2*np.pi) # Displacement of each particle
a = np.zeros_like(i) # 2nd ODE [(d**2)*x_i]/[d*(t**2)], or acceleration, of each particle
v = np.zeros_like(i) # 1st ODE (d*x_i)/(d*t) = u_i, or velocity, of each particle

data = [] # Array to store values of displacements of particles

# Graphic plotting organizing
fig, ax = plt.subplots() # Acceleration against displacement of particles
line, = ax.plot(i, x)


# Calculate & Store the acceleration of all particles exc. two fixed end points
def model(x):
    for j in range(1, n_b - 1):
        # Perform Euler's method w/ the governing equation
        a[j] = (
            (k_young/(rho*h**2)) * 
            (x[j - 1] + x[j + 1] - 2*x[j]) * 
            (1 + alpha*(x[j + 1] - x[j - 1]))
            )
    return a


def graph_from_data(data):
    # Generate animated plotting
    # Function "animate_from_data()" sets particle displacements as y-axis values
    ani = animation.FuncAnimation(
        fig, animate_from_data, frames = int(frame_count), 
        fargs = [data], interval = 0.01, blit = True, repeat = False
        )
    writer = PillowWriter(fps = 25)
    # Save the generated animated plotting
    ani.save("demo_sine{}{}_RK.gif".format(dt,n_b), writer=writer)
    plt.show() # Show the generated animated plotting


# Set particle displacements as y-axis values in plotting
def animate_from_data(frame, data):
    y = data[frame]
    line.set_ydata(y)
    return line,


# Function to write particle displacements at every time frame into a text file
def write_file(data):
    file_name = "FPUT_experiment{}{}_RK.npy".format(dt, n_b)
    # Check whether a file w/ a same name under the file exist or not
    # This ensures particle displacements can be stored cleanly
    if os.path.exists(file_name): # If a same-name exists, remove the previous file
        os.remove(file_name) # Perform clean over-writing
    np.save(file = file_name, arr = data) # Save particle displacement values for plotting


# Program mainframe starts here
# Perform Runge-Kutta method of coefficient computation @ 4th order
for i in range(frame_count):
    data.append(x)

    # Compute the four step slopes
    y_dummy = x
    k1 = model(y_dummy)

    y_dummy = x + k1 * dt/2
    k2 = model(y_dummy)

    y_dummy = x + k2 * dt/2
    k3 = model(y_dummy)

    y_dummy = x + k3 * dt
    k4 = model(y_dummy)

    # Compute the weighted avaerage slope from the four step slopes
    slope = 1/6*k1 + 2/6*k2 + 2/6*k3 + 1/6*k4

    # Compute the velocity & the displacement of particles
    v += dt*slope
    x += dt*v

write_file(data)
print('file saved')
graph_from_data(data)