# DEV: M-H-Fahrudin, S-Li, L-Parker, B-Beale, and P-Wakely-Skinnarl
# Build: v2.0
# Usage: FPUT problem solving program

#------------------------------------------------------------------------------------------------------------

# STARTS: Blackbox input parameters user interface

# Welcome
# To avoid errors, please input acceptable parameter values
# Save this file before run any linked further programs

# String length: L = <POSITIVE INTEGER>
L = 1

# Total simulation time: final_time = <POSITIVE INTEGER>
# Unit: Second
final_time = 12000

# Time step: dt = <POSITIVE FLOAT>
# Unit: Second
dt = 0.1

# Runge-Kutta K steps: RK_n = <POSITIVE INTEGER>
RK_n = 4

# Number of oscillators: n_b = <POSITIVE INTEGER>
n_b = 65

# Young’s modulus: k_young = <POSITIVE FLOAT>
k_young = 0.1

# Mass per unit length of the string: rho = <POSITIVE FLOAT>
rho = 400

# Non-linear coefficient: alpha = <POSITIVE FLOAT>
alpha = 0.25

# Amplitude of initial condition: pert_ini = <POSITIVE FLOAT>
pert_ini = 1.0

# This variable is not required
# Initial condition setting:-
    # "init_type = 1": single sine
    # "init_type = 2": half sine
    # "init_type = 3": parabola
init_type = 3

# Method of solving the Second Order ODE:-
    # "method_type = 1" : Euler Method
    # "method_type = 2" : 4th Order Runge-Kutta Method
method_type = 1

# Do you want to run the simulation with tolerance? (1=YES / 2=NO)
toggle_tol = 2

# The percentage error for tolerance calculation: err = <POSITIVE FLOAT>
# (Only available if simulation is run with tolerance)
err = 1



# ENDS: Blackbox input parameters user interface


#------------------------------------------------------------------------------------------------------------


# STARTS: Whitebox input parameters user interface
    # Call the linked python program
    # Automatically pop out generated plotting

def whitebox_interface():
    
    # Use minor time delays to avoid prompts appear in incorrect sequences
    global time
    import time
    
    print("")
    print("DEV: M-H-Fahrudin, S-Li, L-Parker, B-Beale, P-Wakley-Skinnarland")
    print("Build: v2.0")
    print("Usage: FPUT problem solving program")
    print(""); print("-----"); print("")
    
    print("Welcome"); print("")
    
    print("To run this program in whitebox mode, input 1 to proceed;")
    time.sleep(0.2)
    initiation = input(str("Or input any other character to terminate: "))
    print(""); print("-----")
    
    if initiation == "1":
        
        print(""); print("### PROGRAM INITIATED. ###")
        
        prompt_parameter_array = [
            "String Length", 
            "Total Simulation Time", 
            "Time Step", 
            "Runge-Kutta Steps", 
            "Number of Oscillators", 
            "Young's Modulus", 
            "Mass per Unit Length of the String", 
            "Non-Linear Coefficient", 
            "Amplitude of Initial Condition"
            ]
        
        data_type_parameter_list = [
            0.1, 1, 0.1, 1, 1, 0.1, 0.1, 0.1, 0.1
            ]
        
        global parameter_list
        parameter_list = []
        
        for i in range(len(prompt_parameter_array)):
            print("")
            
            parameter_list.append(
                exceptional_handling(
                    prompt_parameter_array[i], data_type_parameter_list[i], i)
                                  )
      
        print(""); print("-----"); print("")
        print("Please final check the input parameters:"); print("")
        
        for k in range(len(prompt_parameter_array)):
            print(
                prompt_parameter_array[k], ":", parameter_list[k]
                )
        
        print(""); print("-----"); print("")
        
        # Summon the solution program
        # Insert the following program's name here, w/o the suffix ".py"
        
        print(""); print("### PROGRAM EXECUTED. ###")
        
    else:
        print(""); print("### PROGRAM ELIMINATED. ###")

# Exceptional haldling for incorrect value inputs
def exceptional_handling(prompt, data_type, j):
    
    while True:
        try:
            print("-----"); print("")
            print("Please input parameter:", prompt)
            if j == (1 or 2):
                print("Unit: Second")
            
            if type(data_type) is int:
                print(prompt, "= <POSITIVE INTEGER>")
                print("")
                time.sleep(0.2)
                parameter = int(input("Type in here: "))
            elif type(data_type) is float:
                print(prompt, "= <POSITIVE FLOAT>")
                print("")
                time.sleep(0.2)
                parameter = float(input("Type in here: "))
            
            if parameter <=0:
                parameter = int("Error Trigger") # Non-positive input handling
            break
        except ValueError:
            print("### INPUT VALUE ERROR. ###"); print("")
            
    return parameter

whitebox_interface()

# ENDS: Whitebox input parameters user interface