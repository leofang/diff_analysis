import scipy.io
from alignment import get_edge_pos_1D

mat = scipy.io.loadmat("../test_data/unaligned_data.mat")
result1 = get_edge_pos_1D(mat)
result2 = get_edge_pos_1D(mat, fit=True)
print(result1)
print(result2)
