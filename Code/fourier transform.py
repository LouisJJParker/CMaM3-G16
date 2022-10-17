import numpy as np
import matplotlib.pyplot as plt
import start_conditions
import scipy

file_name = "FPUT_experiment0.165.npy"
data = np.load(file_name)
i = np.linspace(0,start_conditions.L,start_conditions.n_b)

def prepare_data(data):
    data