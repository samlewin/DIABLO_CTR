import numpy as np
import h5py 

def create_uniform_grid(NY, domain_size, savedir=''):
    spacing = domain_size/(NY-1)
    y = np.linspace(-domain_size/2-spacing/2, domain_size/2+spacing/2, NY+1)
    grid = h5py.File(savedir +'grid.h5','w')
    grid.create_group('grids')
    grid.create_dataset('grids/y', data=y)
    grid.close()
    if len(savedir) == 0:
        print('Saved file grid.h5 in current working directory')
    else:
        print('Saved file grid.h5 in ' + savedir)
    grid_gyf = (y[1:] + y[:-1])/2
    np.save('grid_gyf', grid_gyf) 
    return grid_gyf

def calc_NY_MPI(NY, NPROCS_Y):
    return (NY-1)/NPROCS_Y+1