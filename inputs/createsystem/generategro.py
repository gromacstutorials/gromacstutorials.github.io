#!/usr/bin/env python
# coding: utf-8

import numpy as np
from molecules import SO4_ion, Na_ion, H20_molecule
from utils import generate_random_location, search_closest_neighbor

# define the box size
Lx, Ly, Lz = [3.36]*3
box = np.array([Lx, Ly, Lz])

# calculate the number of ions
Mh2o = 0.018053 # kg/mol - water
ntotal = 720 # total number of molecule
c = 1.5 # desired concentration in mol/L
nion = c*ntotal*Mh2o/(3*(1+Mh2o*c)) # desired number for the SO4 ion
nwater = ntotal - 3*nion # desired number of water

# define cutoff 
dSO4 = 0.45
dNa = 0.28
dSol = 0.28

# initialized
cpt_residue = 0
cpt_atoms = 0
cpt_SO4 = 0
cpt_Na = 0
cpt_Sol = 0
all_positions = []
all_resnum = []
all_resname = []
all_atname = []
all_attype = []

# add SO4 randomly
atpositions, attypes, resname, atnames = SO4_ion()
while cpt_SO4 < np.int32(nion):
    x_com, y_com, z_com = generate_random_location(box)
    d = search_closest_neighbor(np.array(all_positions), atpositions + np.array([x_com, y_com, z_com]), box)
    if d < dSO4:
        add_residue = False
    else:
        add_residue = True
    if add_residue == True:
        cpt_SO4 += 1
        cpt_residue += 1
        for atposition, attype, atname in zip(atpositions, attypes, atnames):
            cpt_atoms += 1
            x_at, y_at, z_at = atposition
            all_positions.append([x_com+x_at, y_com+y_at, z_com+z_at])
            all_resnum.append(cpt_residue)
            all_resname.append(resname)
            all_atname.append(atname)
            all_attype.append(attype)

# add Na randomly
atpositions, attypes, resname, atnames = Na_ion()
while cpt_Na < np.int32(nion*2):
    x_com, y_com, z_com = generate_random_location(box)
    d = search_closest_neighbor(np.array(all_positions), atpositions + np.array([x_com, y_com, z_com]), box)
    if d < dNa:
        add_residue = False
    else:
        add_residue = True
    if add_residue == True:
        cpt_Na += 1
        cpt_residue += 1
        for atposition, attype, atname in zip(atpositions, attypes, atnames):
            cpt_atoms += 1
            x_at, y_at, z_at = atposition
            all_positions.append([x_com+x_at, y_com+y_at, z_com+z_at])
            all_resnum.append(cpt_residue)
            all_resname.append(resname)
            all_atname.append(atname)
            all_attype.append(attype)

# add water randomly
atpositions, attypes, resname, atnames = H20_molecule()
for x_com in np.arange(dSol/2, Lx, dSol):
    for y_com in np.arange(dSol/2, Ly, dSol):
        for z_com in np.arange(dSol/2, Lz, dSol):
            d = search_closest_neighbor(np.array(all_positions), atpositions + np.array([x_com, y_com, z_com]), box)
            if d < dSol:
                add_residue = False
            else:
                add_residue = True
            if (add_residue == True) & (cpt_Sol < np.int32(nwater)):
                cpt_Sol += 1
                cpt_residue += 1
                for atposition, attype, atname in zip(atpositions, attypes, atnames):
                    cpt_atoms += 1
                    x_at, y_at, z_at = atposition
                    all_positions.append([x_com+x_at, y_com+y_at, z_com+z_at])
                    all_resnum.append(cpt_residue)
                    all_resname.append(resname)
                    all_atname.append(atname)
                    all_attype.append(attype)
            if cpt_Sol >= np.int32(nwater):
                break
print(cpt_Sol, 'out of', np.int32(nwater), 'water molecules created')

# print a few information
print('Lx = '+str(Lx)+' nm, Ly = '+str(Ly)+' nm, Lz = '+str(Lz)+' nm')
print(str(cpt_Na)+' Na ions') 
print(str(cpt_SO4)+' SO4 ions')
print(str(cpt_Sol)+' Sol mols')
Vwater = cpt_Sol/6.022e23*0.018 # kg or litter
Naddion = (cpt_Na+cpt_SO4)/6.022e23 # mol
cion = Naddion/Vwater
print('The ion concentration is '+str(np.round(cion,2))+' mol per litter')

# write conf.gro
f = open('conf.gro', 'w')
f.write('Na2SO4 solution\n')
f.write(str(cpt_atoms)+'\n')
cpt = 0
for resnum, resname, attype, position  in zip(all_resnum, all_resname, all_attype, all_positions):
    x, y, z = position
    cpt += 1
    f.write("{: >5}".format(str(resnum))) # residue number (5 positions, integer) 
    f.write("{: >5}".format(resname)) # residue name (5 characters) 
    f.write("{: >5}".format(attype)) # atom name (5 characters) 
    f.write("{: >5}".format(str(cpt))) # atom number (5 positions, integer)
    f.write("{: >8}".format(str("{:.3f}".format(x)))) # position (in nm, x y z in 3 columns, each 8 positions with 3 decimal places)
    f.write("{: >8}".format(str("{:.3f}".format(y)))) # position (in nm, x y z in 3 columns, each 8 positions with 3 decimal places) 
    f.write("{: >8}".format(str("{:.3f}".format(z)))) # position (in nm, x y z in 3 columns, each 8 positions with 3 decimal places) 
    f.write("\n")
f.write("{: >10}".format(str("{:.5f}".format(Lx)))) # box size
f.write("{: >10}".format(str("{:.5f}".format(Ly)))) # box size
f.write("{: >10}".format(str("{:.5f}".format(Lz)))) # box size
f.write("\n")
f.close()

# write topol.top
f = open('topol.top', 'w')
f.write('#include "ff/forcefield.itp"\n')
f.write('#include "ff/h2o.itp"\n')
f.write('#include "ff/na.itp"\n')
f.write('#include "ff/so4.itp"\n\n')
f.write('[ System ]\n')
f.write('Na2SO4 solution\n\n')
f.write('[ Molecules ]\n')
f.write('SO4 '+ str(cpt_SO4)+'\n')
f.write('Na '+ str(cpt_Na)+'\n')
f.write('SOL '+ str(cpt_Sol)+'\n')
f.close()

