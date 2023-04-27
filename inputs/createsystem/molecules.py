import numpy as np

# define SO4 ion
def SO4_ion():
    Position = np.array([[1.238,   0.587,   1.119], \
        [0.778,   1.501,  -1.263], \
        [-0.962,   1.866,   0.623], \
        [-0.592,  -0.506,  -0.358],\
        [0.115,   0.862,   0.030]])
    Type = ['O1', 'O2', 'O3', 'O4', 'S']
    Resname = 'SO4'
    return Position, Type, Resname

# define Na ion
def Na_ion():
    Position = np.array([[0, 0, 0]])
    Type = ['Na']
    Resname = 'Na'
    return Position, Type, Resname

def H20_molecule():
    Position = np.array([[ 0.    ,  0.    ,  0.    ], \
        [ 0.5858,  0.757 ,  0.    ], \
        [ 0.5858, -0.757 ,  0.    ], \
        [ 0.104 ,  0.    ,  0.    ]])
    Type = ['OW', 'HW1', 'HW2', 'MW']
    Resname = 'Sol'
    return Position, Type, Resname