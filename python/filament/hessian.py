#import CubicEquationSolver
import numpy as np
from numpy.linalg import eig

def hessian_element(field, i, j, dds):
        #shape_extend = np.array(field.shape)+2
        #field_copy = np.zeros(shape_extend)
        dim = field.ndim
        dims = field.shape
        slcLeft = slice(None,-2)
        slcRight = slice(2,None)
        slcMid = slice(1,-1)
        slcMid3D = tuple([slcMid]*dim)
        #field_copy[slcMid3D] = field
        #field_copy[0,:,:]=field_copy[-2,:,:]
        #field_copy[-1,:,:]=field_copy[1,:,:]
        if i==j:
            slcLeft3D = [slcMid]*dim
            slcRight3D = [slcMid]*dim
            slcLeft3D[i] = slcLeft
            slcRight3D[i] = slcRight
            return (field[tuple(slcLeft3D)]+field[tuple(slcRight3D)]-2*field[slcMid3D]) / (dds[i] * dds[i])
        else:
            slcDiag1 = [slcMid]*dim
            slcDiag2 = [slcMid]*dim
            slcDiag3 = [slcMid]*dim
            slcDiag4 = [slcMid]*dim
            slcDiag1[i] = slcRight
            slcDiag1[j] = slcRight
            slcDiag2[i] = slcLeft
            slcDiag2[j] = slcLeft
            slcDiag3[i] = slcRight
            slcDiag3[j] = slcLeft
            slcDiag4[i] = slcLeft
            slcDiag4[j] = slcRight
            return (field[tuple(slcDiag1)] + field[tuple(slcDiag2)] - field[tuple(slcDiag3)] - field[tuple(slcDiag4)]) / (4 * dds[i] * dds[j])

def hessian(array,dds):
       dim = array.ndim
       sh_orig = array.shape
       sh_arr = np.asarray(sh_orig)
       sh_arr = sh_arr-2
       sh = tuple(sh_arr)+(dim,dim)
       res = np.zeros(sh)
       for i in range(dim):
              for j in range(dim):
                        res[...,i,j] = hessian_element(array, i, j, dds)
       return res

def eigensystem(hessian_field):
       return eig(hessian_field)

def sort_hessian(E,V):
       # this finds the ordering on eigenvalues
       ordering_E = E.argsort(axis=-1)
       # applies ordering to eigenvalues
       E[:] = np.take_along_axis(E, ordering_E, axis=-1)
       # extends the ordering (adds dummy axis to the second to last place)
       # needs to add it 3 times (repeat). In general, it needs to be added n times
       # n = E.shape[-1] - this gives the number of eigenvalues, but also the number of components of each eigenvector.
       ordering_V = np.repeat(ordering_E[...,np.newaxis,:], E.shape[-1], axis=-2)
       # applies ordering to V and transposes the last two levels, so that the components span the last index as rows
       V[:] = np.moveaxis(np.take_along_axis(V, ordering_V, axis=-1),-2,-1)

def extend_periodic(field, ghostCount):
       field_dim = field.ndim
       shape_orig = field.shape
       shape_extended = np.array(field.shape) + 2 * ghostCount
       field_copy = np.zeros(shape_extended)
       # coords is a list of individual coordinates of the original field. In a way, it mimicks the whole grid but in indices.
       # it can be used directly with the field at a single position: field[coords[0], coords[1], ...]
       # [[-1, 0, 1, ..., dimX - 1, dimX], [-1, 0, ...], ...]
       # two things (specific for ghostCount = 1): -1 and dimX are added on both ends
       # -1 and dimX are then coverted to their periodic counterpart by adding modulo dimX
       coords = [np.arange(-ghostCount, shape_orig[dim] + ghostCount) % shape_orig[dim] for dim in range(field_dim)]
       # the following needs explanation.
       # asterisk (*) before a thing unpacks the elements of that thing into separate arguments. handy.
       # np.ix_(stuff) converts stuff into a meshgrid tuple, exactly how you'd imagine it.
       # a[np.ix_([1,3],[2,5])] returns the array [[a[1,2] a[1,5]], [a[3,2] a[3,5]]]
       indexing_grid = np.ix_(*coords)
       # as promised, the ix_ object can be stuffed right into the square brackets to span the field.
       # the extraneous -1, dim padded at the end were turned (by modulo) into something that is still in range
       field_copy = field[indexing_grid]
       return field_copy

# goals:
# 2. statistics
#    a) blobs: e1, e2, e3 < 0
#    b) lines: e1, e2 < 0 and |e1|, |e2| > |e3|
#    c) sheets: e1 < 0 and |e1| > |e2|, |e3|
# 3. for lines, do B dot v3 (the one eigenvector that corresponds to the largest eigenvalue)

# (512, 512, 512)
# [[511, 0, 1, ..., 511, 0] mod 512 x [0, 1, ..., 512] x [0, 1, ..., 512]]