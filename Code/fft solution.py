from scipy import fftpack
from scipy.sgnal import find_peaks
import numpy as np
import matplotlib.pyplot as plt

# f = data
# N = number of oscillations
# dt = time step

y = fftpack.rfft(f)
yf = np.abs(y)
yf_amp = yf/(N/2)
yf_amp_1 = yf_amp[:,2]
yf_amp_2 = yf_amp[:,4]
yf_amp_3 = yf_amp[:,6]
yf_amp_4 = yf_amp[:,8]

p1, _ = find_peaks(yf_amp_1)
p1_val = np.take(yf_amp_1, p1)

p2, _ = find_peaks(yf_amp_2)
p2_val = np.take(yf_amp_2, p2)

p3, _ = find_peaks(yf_amp_3)
p3_val = np.take(yf_amp_3, p3)

p4, _ = find_peaks(yf_amp_4)
p4_val = np.take(yf_amp_4, p4)

p1 = p1*dt
p2 = p2*dt
p3 = p3*dt
p4 = p4*dt

plt.xlabel('Time')
plt.ylabel('Energy')
plt.title('FFT')

plt.plot(p1, p1_val, 'r-', label='First Coefficient')
plt.plot(p2, p2_val, 'b-', label='Second Coefficient')
plt.plot(p3, p3_val, 'g-', label='Third Coefficient')
plt.plot(p4, p4_val, 'y-', label='Fourth Coefficient')
plt.legend()
plt.show()
