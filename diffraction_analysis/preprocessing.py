import numpy as np


def canonicalize_1D(inarr, copy=True):
    '''
    Bring the input array of shape (det_row, det_col, num_pos, num_angle)
    to the canonical form (num_angle, num_pos, det_row, det_col).

    Arguments:
        inarr: input array with inarr.ndim==4.
        copy: if True (default), a copy is returned, else a view whenever possible. 

    Return:
        out: a view (or a copy) of the input array in the canonical form.
    '''
    if inarr.ndim != 4:
        raise TypeError("The input should be 4D (i.e. from 1D scan data).")

    det_row, det_col, num_pos, num_angle = inarr.shape
    ## the transpose action is identical to the double for loop below:
    #for i in range(num_angle):
    #    for j in range(num_pos):
    #        out[i, j] = inarr[..., j, i]
    out = np.transpose(inarr, axes=(3, 2, 0, 1))

    if copy:
        return out.copy()
    else:
        return out

def canonicalize_2D(inarr, copy=False):
    '''
    Bring the input array of shape (det_row, det_col, num_posx, num_posy, num_angle)
    to the canonical form (num_angle, num_posx, num_posy, det_row, det_col).

    Arguments:
        inarr: input array with inarr.ndim==5.
        copy: if False (default), a view is returned whenever possible, else a copy.

    Return:
        out: a view (or a copy) of the input array in the canonical form.
    '''
    if inarr.ndim != 5:
        raise TypeError("The input should be 5D (i.e. from 2D scan data).")

    det_row, det_col, num_posx, num_posy, num_angle = inarr.shape
    out = np.transpose(inarr, axes=(4, 2, 3, 0, 1))

    if copy:
        return out.copy()
    else:
        return out
