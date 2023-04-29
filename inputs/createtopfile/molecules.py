import numpy as np

# define SO4 ion
def SO4_ion():
    Position = np.array([[0.1238,  0.0587,   0.1119], \
        [0.0778,   0.1501,  -0.1263], \
        [-0.0962,   0.1866,   0.0623], \
        [-0.0592,  -0.0506,  -0.0358],\
        [0.0115,   0.0862,   0.0030]])
    Type = ['OS', 'OS', 'OS', 'OS', 'SO']
    Name = ['O1', 'O2', 'O3', 'O4', 'S1']
    Resname = 'SO4'
    return Position, Type, Resname, Name

# define Na ion
def Na_ion():
    Position = np.array([[0, 0, 0]])
    Type = ['Na']
    Name = ['Na1']
    Resname = 'Na'
    return Position, Type, Resname, Name

# define water molecule
def H20_molecule():
    Position = np.array([[ 0.    ,  0.    ,  0.    ], \
        [ 0.05858,  0.0757 ,  0.    ], \
        [ 0.05858, -0.0757 ,  0.    ], \
        [ 0.0104 ,  0.    ,  0.    ]])
    Type = ['OW', 'HW', 'HW', 'MW']
    Name = ['OW1', 'HW1', 'HW2', 'MW1']
    Resname = 'Sol'
    return Position, Type, Resname, Name