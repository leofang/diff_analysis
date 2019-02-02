import scipy.io
import numpy as np

from alignment import get_edge_pos_1D
from preprocessing import canonicalize_1D

mat = scipy.io.loadmat("../test_data/unaligned_data.mat")
pos = get_edge_pos_1D(mat) # faster
#pos = get_edge_pos_1D(mat, fit=True)

diff_data = canonicalize_1D(mat['assembled_det_data'])
print(diff_data.shape)

shift = np.int(np.median(pos))-pos
print(shift_max, shift_min)
for i in range(diff_data.shape[0]):
    diff_data[i] = np.roll(diff_data[i], shift[i], axis=0)
shift_max = np.max(shift) # >0
shift_min = np.min(shift) # <0
new_diff_data = diff_data[:, shift_max:shift_min, :, :]
print(new_diff_data.shape)
