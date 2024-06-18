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
    return None

print('Enter NY')
NY = input()
print('Enter domain size')
domain_size = input()

create_uniform_grid(int(NY), int(domain_size))


   