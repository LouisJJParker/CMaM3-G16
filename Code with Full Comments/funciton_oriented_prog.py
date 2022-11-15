# Import essential built-in modules
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import PillowWriter
from scipy.integrate import solve_ivp
import start_conditions
import os


# Import initial conditions from file "start_condition.py"
L = start_conditions.L # String length
n_b = start_conditions.n_b + 2 # Number of oscillators + the two fixed end points
rho = start_conditions.rho # Mass per unit length (Density)
alpha = start_conditions.alpha # Non-linear coefficient
k_young = start_conditions.k_young # Young's modulus
init_type = start_conditions.init_type # Initial conditions

dt = start_conditions.dt # Time step
final_time = start_conditions.final_time # Total simulation time in seconds

frame_count = int(final_time//dt) # Number of simulated frames

h = L/(n_b - 3) # Latice spacing distance between two neighboring particles
i = np.linspace(0, L, n_b) # Indexes of particles
x = np.sin(i*2*np.pi) # Variable of displacement of each particle
a = np.zeros_like(i) # Array storing [(d**2)*x_i]/[d*(t**2)], or acceleration, of each particle
v = np.zeros_like(i) # Array storing (d*x_i)/(d*t) = u_i, or velocity, of each particle

# Graphic plotting organizing
fig, ax = plt.subplots() # Acceleration against displacement of particles
line, = ax.plot(i,x)


# Program main procedure
def main():
    print("FPUT_experiment{}{}.npy".format(dt,n_b))
    # Save the array of particle displacements at every time fram into a text file
    write_file(run_without_graph())
    # Extract file-stored particle displacements to generate plotting
    graph_from_file("FPUT_experiment{}{}.npy".format(dt,n_b))
    #plot_at_t(10000, read_file("FPUT_experiment0.165.npy"))


# Function of extracting displacement, 
def read_file(file_name):
    data = np.load(file = file_name)
    return data


# Function to write particle displacements at every time frame into a text file
def write_file(data):
    file_name = "FPUT_experiment{}{}.npy".format(dt, n_b)
    # Check whether a file w/ a same name under the file exist or not
    # This ensures particle displacements can be stored cleanly
    if os.path.exists(file_name): # If a same-name exists, remove the previous file
        os.remove(file_name) # Perform clean over-writing
    np.save(file = file_name, arr = data) # Save particle displacement values for plotting


# Calculate & Store the acceleration of all particles exc. two fixed end points
# REDUNDANT intake variables: v
def simulate(a, v, x):
    # Perform Euler's method w/ the governing equation
    for j in range(1, n_b - 1):
        # Rearrange m & Recall m = rho*h
        a[j] = ((k_young/(rho*h**2)) * 
                (x[j - 1]+x[j + 1] - 2*x[j]) * 
                (1 + alpha*(x[j + 1] - x[j - 1])))
    return a, v, x


# Use Euler's method to obtain displacements of particles
def euler_method(a,v,x):
    a, v, x = simulate(a,v,x)
    v += a*dt # Obtain the particle velocity at time intervals from acceleration
    x += v*dt # Obtain the displacement at time intervals from velocity
    
    # Return particle displacements, "x", to be stored in array "data[]"
    # Return value "x" will be used in procedure "run_without_graph()"
    return x


# Store the displacement of each numerically at every time frame before graph plotting
def run_without_graph():
    data = [] # Array to store numerical displacement of particles at each time frame
    for i in range(frame_count):
        data.append(euler_method(a,v,x).copy())
    
    # Return to write particle displacements at every time frame into a text file
    return data


# !!! REDUNDANT FUNCTION !!
# REDUNDANT intake variables: t
def animate(t):
    y = euler_method(a,v,x)
    line.set_ydata(y)
    return line,

# Set particle displacements as y-axis values in plotting
def animate_from_data(frame, data):
    y = data[frame]
    line.set_ydata(y)
    return line,


# !!! REDUNDANT PROCEDURE !!!
def show_graph():
    ani = animation.FuncAnimation(
        fig, animate, interval = 0.001, blit = True, save_count = 50)

    plt.show()


# !!! REDUNDANT FUNCTION !!!
def graph_from_data(data):
    ani = animation.FuncAnimation(
        fig, animate_from_data, frames=range(0,frame_count,40), 
        fargs=[data], interval=0.1, blit=True, repeat = False
        )

    plt.show()


# Extract stored particle displacements from the text file to generate plotting
def graph_from_file(file_name):
    data = read_file(file_name) # Extract the file containing particle displacements
    # Generate animated plotting
    # Function "animate_from_data()" sets particle displacements as y-axis values
    ani = animation.FuncAnimation(
        fig, animate_from_data, frames = range(0,frame_count,1), 
        fargs = [data], interval = 0.1, blit = True, repeat = False
        )
    
    writer = PillowWriter(fps = 25)
    # Save the generated animated plotting
    ani.save("demo_sine{}{}.gif".format(dt,n_b), writer = writer)
    plt.show() # Show the generated animated plotting


# !!! REDUNDANT FUNCTION !!!
def plot_at_t(t,data):
    line.set_ydata(data[t])

    plt.show()

# Program initation calling
main()
