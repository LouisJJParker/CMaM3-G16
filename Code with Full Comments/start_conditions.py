# DEV: M-H-Fahrudin, S-Li, L-Parker, B-Beale, and P-Wakely-Skinnarl
# Build: v2.0
# Usage: FPUT problem solving program

# -----

# STARTS: Blackbox input parameters user interface

# Welcome
# To avoid errors, please input acceptable parameter values
# Save this file before run any linked further programs

# String length: L = <POSITIVE FLOAT>
L = 1.0

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

# Initial condition setting
    # "init_type = 1": single sine
    # "init_type = 2": half sine
    # "init_type = 3": parabola
init_type = 1

# Method chosen to solve the Second Order ODE
    # "method_type = 1": Euler's Method
    # "method_type = 2": Runge-Kutta Method on 4th Order
method_type = 1

# Toggle of running w/ tolerance
    # "toggle_tol = 1": YES, w/ tolerance
    # "toggle_tol = 2": NO, w/o tolerance
toggle_tol = 1

# Tolerance percentage error: err = <POSITIVE FLOAT>
    # Only available w/ toerance toggled on
err = 10

# ENDS: Blackbox input parameters user interface


# -----


# STARTS: Whitebox input parameters user interface
    # Call the linked python program
    # Automatically pop out generated plotting

def whitebox_interface():
    
    # Use minor time delays to avoid prompts appear in incorrect sequences
    global time
    import time
    
    print(
        "n\DEV: M-H-Fahrudin, S-Li, L-Parker, B-Beale, P-Wakely-Skinnarland")
    print("Build: v2.0")
    print("Usage: FPUT problem solving program\n")
    print("-----\n")
    
    print("Welcome\n")
    
    print("To run this program in whitebox mode, input 1 to proceed:")
    time.sleep(0.2)
    initiation = input(str("Or input any other character to terminate: "))
    print("\n-----")
    
    if initiation == "1":
        
        print("\n### PROGRAM INITIATED. ###")
        
        prompt_parameter_array = [
            "String Length", 
            "Total Simulation Time", 
            "Time Step", 
            "Runge-Kutta Steps", 
            "Number of Oscillators", 
            "Young’s Modulus", 
            "Mass per Unit Length of the String", 
            "Non-Linear Coefficient", 
            "Amplitude of Initial Condition",
            "Initial Condition Setting",
            "Solution Method Setting",
            "Tolerance ON/OFF Switch",
            "Tolerance Percentage Error"
            ]
        
        data_type_parameter_list = [
            0.1, 1, 0.1, 1, 1, 0.1, 0.1, 0.1, 0.1, 1, 1, 1, 0.1
            ]
        
        parameter_list = []
        
        i = 0
        tolerance_ON = True
        
        while (i <= (len(prompt_parameter_array) - 1)) and (
                tolerance_ON == True
                ):
            print("")
            
            parameter_list.append(
                exceptional_handling(
                    prompt_parameter_array[i], data_type_parameter_list[i], i
                    ))
            
            if (i == len(prompt_parameter_array) - 2) and (
                    parameter_list[i] == "NO, continue w/o tolerance"):
                tolerance_ON = False
            else:
                i += 1
      
        print("\n-----\n")
        print("### Parameters inserted: ###\n")
        
        for k in range(len(parameter_list)):
            print(
                prompt_parameter_array[k], ":", parameter_list[k]
                )
        
        print("\n-----\n")
        
        # Parameter extraction before calling the solution prgram
        L = parameter_list[0]
        final_time = parameter_list[1]
        dt = parameter_list[2]
        RK_n = parameter_list[3]
        n_b = parameter_list[4]
        k_young = parameter_list[5]
        rho = parameter_list[6]
        alpha = parameter_list[7]
        pert_ini = parameter_list[8]
        
        if parameter_list[9] == "Single-sine initial condition":
            init_type = 1
        elif parameter_list[9] == "Half-sine initial condition":
            init_type = 2
        else:
            init_type = 3
        
        if parameter_list[10] == "Euler's Method":
            method_type = 1
        else:
            method_type = 2
        
        if parameter_list[11] == "YES, continue w/ tolerance":
            toggle_tol = 1
        else:
            toggle_tol = 2
        
        if toggle_tol == 1:
            err = parameter_list[12]
        
        
        # Summon the solution program
        # Insert the following program's name here, w/o the suffix ".py"
        import program_to_be_called
        program_to_be_called
        
        print("\n### PROGRAM EXECUTED. ###")
        
    else:
        print("\n### PROGRAM ELIMINATED. ###")

# Exceptional haldling for incorrect value inputs
def exceptional_handling(prompt, data_type, j):
    
    while True:
        try:
            print("-----\n")
            print("Please input parameter:", prompt)
            if j == (1 or 2):
                print("Unit: Second")
            
            # Special case when choosing the Initial Condition
            if j == 9:
                init_condition_prompt = [
                    "Single-sine initial condition",
                    "Half-sine initial condition",
                    "Parabolic initial condition"
                    ]
                print(prompt, "= Either of {1 OR 2 OR 3}\n")
                
                for k in range(len(init_condition_prompt)):
                    print((k + 1), ":", init_condition_prompt[k])

                print("")
                time.sleep(0.2)
                parameter = int(input("Type in here: "))
                
                # Deny out-of-range option input
                if (parameter not in [1, 2, 3]) == True:
                    parameter = int("Error trigger") # Re-input trigger
            
            # Special case when choosing the Solution Method Type
            elif j == 10:
                solution_method_prompt = [
                    "Euler's Method", "Runge-Kutta Method on 4th Order"
                    ]
                print(prompt, "= Either of {1 OR 2}\n")
                
                for k in range(len(solution_method_prompt)):
                    print((k + 1), ":", solution_method_prompt[k])
                
                print("")
                time.sleep(0.2)
                parameter = int(input("Type in here: "))
                
                # Deny out-of-range option input
                if (parameter not in [1, 2]) == True:
                    parameter = int("Error trigger") # Re-input trigger
                
            # Special case when toggling the Tolerance ON/OFF
            elif j == 11:
                tolerance_toggle_prompt = [
                    "YES, continue w/ tolerance",
                    "NO, continue w/o tolerance"
                    ]
                print(prompt, "= Either of {1 OR 2}\n")
                
                for k in range(len(tolerance_toggle_prompt)):
                    print((k + 1), ":", tolerance_toggle_prompt[k])
                
                print("")
                time.sleep(0.2)
                parameter = int(input("Type in here: "))
                
                # Deny out-of-range option input
                if (parameter not in [1, 2]) == True:
                    parameter = int("Error trigger") # Re-input trigger
                
            else:
                if type(data_type) is int:
                    print(prompt, "= <POSITIVE INTEGER>\n")
                    time.sleep(0.2)
                    parameter = int(input("Type in here: "))
                elif type(data_type) is float:
                    print(prompt, "= <POSITIVE FLOAT>\n")
                    time.sleep(0.2)
                    parameter = float(input("Type in here: "))
                
                # Deny non-positive input
                if parameter <=0:
                    parameter = int("Error Trigger") # Re-input trigger
            break
        except ValueError:
            print("### INPUT VALUE ERROR. ###\n")
    
    if j == 9:
        parameter = init_condition_prompt[parameter - 1]
    elif j == 10:
        parameter = solution_method_prompt[parameter - 1]
    elif j == 11:
        parameter = tolerance_toggle_prompt[parameter - 1]
            
    return parameter

whitebox_interface()

# ENDS: Whitebox input parameters user interface
