'''
November 2022

Computational Method and Modelling 3

Project Group 16

Fermi-Pasta-Ulam-Tsingou experiment simulation

Euler Method of solving the problem
'''
#------------------------------------------------------------------------------------------------------------

# Import essential built-in modules
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import PillowWriter
from scipy.integrate import solve_ivp
import start_conditions
import os
from scipy import fftpack
from scipy.signal import find_peaks
import random

#------------------------------------------------------------------------------------------------------------

# Import initial conditions from file "start_condition.py"
L = start_conditions.L # String length
n_b = start_conditions.n_b # Number of oscillators + the two fixed end points
rho = start_conditions.rho # Mass per unit length (Density)
alpha = start_conditions.alpha # Non-linear coefficient
k_young = start_conditions.k_young # Young's modulus
init_type = start_conditions.init_type # Initial profile displacement
pert_ini = start_conditions.pert_ini # Amplitude of initial condition

dt = start_conditions.dt # Time step
final_time = start_conditions.final_time # Total simulation time in seconds

method_type = start_conditions.method_type

toggle_tol = start_conditions.toggle_tol
err = start_conditions.err

#------------------------------------------------------------------------------------------------------------

# Initial calculations before running the simulations
frame_count = math.ceil(final_time//dt) # Number of simulated frames

h = L/(n_b - 3) # Latice spacing distance between two neighboring particles
i = np.linspace(0, L, n_b) # Indexes of particles
x = np.zeros_like(i) # Array storing the displacement of each partice
a = np.zeros_like(i) # Array storing [(d**2)*x_i]/[d*(t**2)], or acceleration, of each particle
v = np.zeros_like(i) # Array storing (d*x_i)/(d*t) = u_i, or velocity, of each particle


# Creating the initial displacement of each particle based on the chosen initial profile displacement
if init_type == 1:
    x = (pert_ini)*np.sin(2*np.pi*i) # Gives out the single sine wave diplacement profile
elif init_type == 2:
    x = (pert_ini)*np.sin(np.pi*i) # Gives out the half sine wave displacement profile
elif init_type == 3:
    x = -(pert_ini)*np.pi(np.pi*i) # Gives out the parabola displacement profile

x[0] = x[-1] = 0

# Graphic plotting organizing
fig, ax = plt.subplots() # Acceleration against displacement of particles
line, = ax.plot(i,x)

#------------------------------------------------------------------------------------------------------------

r = [-1, 1]

def random_number_generator():
    f = random.choice(r)
    fx = f*random.random()
    return fx

tol = np.zeros_like(i)
for j in range(n_b):
    tol[j] = (100 + err * (random_number_generator()))/100

tol[0] = 0

#------------------------------------------------------------------------------------------------------------
'''
Functions created for solving the problems and running the simulations
'''
# Program main procedure
def main():
    print('###  MAIN PROGRAM INITIATED ###')
    print("\nFile Name = FPUT_experiment_{}dt_{}N_ver{}.{}.npy".format(dt,n_b,method_type,toggle_tol))
    
    if toggle_tol == 1:
        print("\n   Tolerance = ON ")
    elif toggle_tol == 2:
        print("\n   Tolerance = OFF ")
    
    if method_type == 1:
        print("   Solving with Euler Method \n")
    elif method_type == 2:
        print("   Solving witha 4th Order Runge-Kutta Method \n")
        
    # Save the array of particle displacements at every time fram into a text file
    print('Initiating calculations \n')
    write_file(run_without_graph())
    print('Calculations finshed. File saved. \n')
    
    # Extract file-stored particle displacements to generate plotting
    print('Initiating simulation video creation \n')
    graph_from_file("FPUT_experiment_{}dt_{}N_ver{}.{}.npy".format(dt,n_b,method_type,toggle_tol))
    print('Video created and saved \n')
    print("File Name = Simulation_{}dt_{}N_ver{}{}.gif\n".format(dt,n_b,method_type,toggle_tol))
    
    print("Initiating plotting at certain time\n")
    plot_at_t("FPUT_experiment_{}dt_{}N_ver{}.{}.npy".format(dt,n_b,method_type,toggle_tol))
    print("Plotting finished and shown\n")
    
    # Plots the Fast Fourier Transform of the simulation
    print('Initiating FFT modal analysis \n')
    FFT("FPUT_experiment_{}dt_{}N_ver{}.{}.npy".format(dt,n_b,method_type,toggle_tol))
    print("FFT solution plotted \n")
    print('### MAIN PROGRAM HAS ENDED ###')

#------------------------------------------------------------------------------------------------------------

# Function of extracting displacement, 
def read_file(file_name):
    data = np.load(file = file_name)
    return data


# Function to write particle displacements at every time frame into a text file
def write_file(data):
    file_name = "FPUT_experiment_{}dt_{}N_ver{}.{}.npy".format(dt, n_b, method_type,toggle_tol)
    # Check whether a file w/ a same name under the file exist or not
    # This ensures particle displacements can be stored cleanly
    if os.path.exists(file_name): # If a same-name exists, remove the previous file
        os.remove(file_name) # Perform clean over-writing
    np.save(file = file_name, arr = data) # Save particle displacement values for plotting

#------------------------------------------------------------------------------------------------------------

# Calculate & Store the acceleration of all particles exc. two fixed end points
def model(x):
    if toggle_tol == 1:
        for j in range(1, n_b - 1):
            # Perform Euler's method w/ the governing equation
            a[j] = (
                (k_young/(rho*h**2)) * 
                (x[j - 1]*tol[j - 1] + x[j + 1]*tol[j + 1] - 2*x[j]*tol[j]) * 
                (1 + alpha*(x[j + 1]*tol[j + 1] - x[j - 1]*tol[j - 1]))
                )
    
    elif toggle_tol == 2:
        for j in range(1, n_b - 1):
            # Perform Euler's method w/ the governing equation
            a[j] = (
                (k_young/(rho*h**2)) * 
                (x[j - 1] + x[j + 1] - 2*x[j]) * 
                (1 + alpha*(x[j + 1] - x[j - 1]))
                )
    return a

#------------------------------------------------------------------------------------------------------------

# Use Euler's method to obtain displacements of particles
def euler_method(a,v,x):
    a = model(x)
    v += a*dt # Obtain the particle velocity at time intervals from acceleration
    x += v*dt # Obtain the displacement at time intervals from velocity
    
    # Return particle displacements, "x", to be stored in array "data[]"
    # Return value "x" will be used in procedure "run_without_graph()"
    return x

#------------------------------------------------------------------------------------------------------------

def RK_method(a,v,x):
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
    a = 1/6*k1 + 2/6*k2 + 2/6*k3 + 1/6*k4

    # Compute the velocity & the displacement of particles
    v += dt*a
    x += dt*v
    
    return x

#------------------------------------------------------------------------------------------------------------

# Store the displacement of each numerically at every time frame before graph plotting
def run_without_graph():
    data = [] # Array to store numerical displacement of particles at each time frame
    if method_type == 1:
        for i in range(frame_count):
            data.append(euler_method(a,v,x).copy())
    
    elif method_type == 2:
        for i in range(frame_count):
            data.append(RK_method(a,v,x).copy())
    # Return to write particle displacements at every time frame into a text file
    return data

# Set particle displacements as y-axis values in plotting
def animate_from_data(frame, data):
    y = data[frame]
    line.set_ydata(y)
    return line,

# Extract stored particle displacements from the text file to generate plotting
def graph_from_file(file_name):
    data = read_file(file_name) # Extract the file containing particle displacements
    # Generate animated plotting
    # Function "animate_from_data()" sets particle displacements as y-axis values
    ani = animation.FuncAnimation(
        fig, animate_from_data, frames = range(0,frame_count,50), 
        fargs = [data], interval = 0.1, blit = True, repeat = False
        )
    
    writer = PillowWriter(fps = 80)
    # Save the generated animated plotting
    ani.save("Simulation_{}dt_{}N_ver{}.{}.gif".format(dt,n_b,method_type,toggle_tol), writer = writer)
    plt.show() # Show the generated animated plotting

#------------------------------------------------------------------------------------------------------------

def plot_at_t(file_name):
    data = read_file(file_name)
    
    print("Asking for the user's input to plot at certain time")
    
    t1 = int(input("Enter the first time plot: "))
    t10 = math.ceil(t1//dt)
    t2 = int(input("Enter the second time plot: "))
    t20 = math.ceil(t2//dt)
    t3 = int(input("Enter the third time plot: "))
    t30 = math.ceil(t3//dt)
    
    print("\nPlotting for all time plots")
    
    plt.title("Plot at time {}".format(t1))
    plt.ylim(-pert_ini, pert_ini)
    plt.plot(i, data[t10])
    plt.ylabel('Displacement of Oscillators')
    plt.show()
    
    plt.title("Plot at time {}".format(t2))
    plt.ylim(-pert_ini, pert_ini)
    plt.plot(i, data[t20])
    plt.ylabel('Displacement of Oscillators')
    plt.show()
    
    plt.title("Plot at time {}".format(t3))
    plt.ylim(-pert_ini, pert_ini)
    plt.plot(i, data[t30])
    plt.ylabel('Displacement of Oscillators')
    plt.show()
       
#------------------------------------------------------------------------------------------------------------ 
 
def FFT(file_name):
    data = read_file(file_name) # Extract the file containing particle displacements
    
    # By using the available fft module 'fftpack' from scipy, we are able to obtain the fourier transform of our simulation
    
    y = fftpack.rfft(data) # finding the fourier transfrom of the data
    yf = np.abs(y) # only taking the absolute real value returned by the fourier done
    yf_amp = yf/(n_b/2) # calculation done to find the amplitude
    yf_amp_1 = yf_amp[:,2] # slicing the array to obtain only the first coefficient amplitude
    yf_amp_2 = yf_amp[:,4] # slicing to obtain the second coefficient amplitude
    yf_amp_3 = yf_amp[:,6] # slicing to obtain the third coefficient amplitude
    yf_amp_4 = yf_amp[:,8] # slicing to obtain the fourth coefficient amplitude
    
    # Using the find_peaks function imported, find the peaks at each slicing

    p1, _ = find_peaks(yf_amp_1)   # Finding the time at which the first peak is
    p1_val = np.take(yf_amp_1, p1) # Taking the value at that time

    p2, _ = find_peaks(yf_amp_2)   # Finding the time at which the second peak is
    p2_val = np.take(yf_amp_2, p2) # Taking the value at that time

    p3, _ = find_peaks(yf_amp_3)   # Finding the time at which the third peak is
    p3_val = np.take(yf_amp_3, p3) # Taking the value at that time

    p4, _ = find_peaks(yf_amp_4)   # Finding the time at which the fourth peak is
    p4_val = np.take(yf_amp_4, p4) # Taking the value at that time
    
    # Changing the time step vaue into time
    p1 = p1*dt
    p2 = p2*dt
    p3 = p3*dt
    p4 = p4*dt
    
    # Labeling the graph
    plt.xlabel('Time')
    plt.ylabel('Energy')
    plt.title('FFT')
    
    # PLotting the all the FFT solution on a single figure
    # Label each plots and placing a legend tab
    plt.plot(p1, p1_val, 'r-', label='First Coefficient')
    plt.plot(p2, p2_val, 'b-', label='Second Coefficient')
    plt.plot(p3, p3_val, 'g-', label='Third Coefficient')
    plt.plot(p4, p4_val, 'y-', label='Fourth Coefficient')
    plt.legend()
    plt.show()

#------------------------------------------------------------------------------------------------------------

# Program initation calling
main()