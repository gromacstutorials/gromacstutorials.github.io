
import numpy as np
from numpy.linalg import norm

def generate_random_location(box):
    """Generate a random location within a given box."""  
    return np.random.rand(3)*box

def search_closest_neighbor(XYZ_neighbor, XYZ_molecule, box):
    """Search neighbor in a box and return the closest distance.
        
    If the neighbor list is empty, then the box size is returned.
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
