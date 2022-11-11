import numpy as np
import matplotlib.pyplot as plt
import start_conditions
from scipy.fft import fft, fftfreq, fftshift, ifft, irfft, rfft, rfftfreq, dst, idst
import matplotlib.animation as animation
from matplotlib.animation import PillowWriter

file_name = "FPUT_experiment0.167.npy"
data = np.load(file_name)
frame_count = len(data)


#need to perform fft to data[t] for t in range(len(data))
#sampling rate = 1/i I think, therefore sr = 1/67


fig, axs = plt.subplots()

# line_1, = axs[0].plot(i,x)
# starting_data = abs(rfft(data[0])/(len(data[0])/2))
# line, = axs.plot(starting_data)
# axs.set_xlim(0,5)
# graph_from_data(data)
# plt.show()

y = []
for i in data:
    y.append(np.abs(rfft(i)/(len(data[0])/2))[:5])

y = np.array(y).T
time = np.linspace(0,start_conditions.final_time, len(y[0]))

for i in range(5):
    axs.plot(time, y[i])

plt.show()