import scipy.io
import numpy as np
import matplotlib.pyplot as plt

from alignment import get_edge_pos_1D
from preprocessing import canonicalize_1D

mat = scipy.io.loadmat("../test_data/unaligned_data.mat")
pos = get_edge_pos_1D(mat) # faster
#pos = get_edge_pos_1D(mat, fit=True)

#flures = mat['fluorescence'].swapaxes(0,1)[:, shift_max:shift_min]
flures = mat['fluorescence'].swapaxes(0,1)
diff_data = canonicalize_1D(mat['assembled_det_data'])
print(diff_data.shape)

shift = np.int(np.median(pos))-pos
for i in range(diff_data.shape[0]):
    diff_data[i] = np.roll(diff_data[i], shift[i], axis=0)
    flures[i] = np.roll(flures[i], shift[i], axis=0)
shift_max = np.max(shift) # >0
shift_min = np.min(shift) # <0
print(shift_max, shift_min)
new_flures = flures[:, shift_max:shift_min]
new_diff_data = diff_data[:, shift_max:shift_min, :, :]
print(new_diff_data.shape)

def comp(i, flures, diff):
    line = flures[i]
    line2 = np.sum(diff[i], axis=(1,2))
    ave = np.mean(line)
    ave2 = np.mean(line2)
    plt.plot(line/ave*ave2)
    plt.plot(line2)

plt.figure()
for i in range(new_flures.shape[0]):
    plt.plot(new_flures[i])

for i in range(0, 31, 5):
    plt.figure()
    comp(i, new_flures, new_diff_data)
    #comp(i, flures, diff_data)
plt.show()
