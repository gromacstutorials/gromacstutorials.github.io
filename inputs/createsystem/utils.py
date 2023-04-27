
import numpy as np
from numpy.linalg import norm

def generate_random_location(Lx, Ly, Lz):
    """
    generate a random location within a given box
    """  
    x = np.random.rand()*Lx
    y = np.random.rand()*Ly
    z = np.random.rand()*Lz
    return x, y, z

def search_closest_neighbor(XYZ_neighbor, XYZ_molecule, box):
    """
        Search neighbor in a box and return the closest distance with a molecule
        Periodic boundary conditions are automatically accounted
    """
    if len(np.array(XYZ_neighbor)) == 0:
        min_distance = np.max(box)
    else:
        min_distance = np.max(box)
        for XYZ_atom in XYZ_molecule:
            dxdydz = np.remainder(XYZ_neighbor - XYZ_atom + box/2., box) - box/2.
            min_distance = np.min([min_distance,np.min(norm(dxdydz,axis=1))])
    return min_distance