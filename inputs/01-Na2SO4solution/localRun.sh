#!/bin/sh

gmx grompp -f input/min.mdp -c conf.gro -p topol.top -o min -pp min -po min
gmx mdrun -v -deffnm min

gmx grompp -f input/nvt.mdp -c min.gro -p topol.top -o nvt -pp nvt -po nvt
gmx mdrun -v -deffnm nvt

gmx grompp -f input/npt.mdp -c nvt.gro -p topol.top -o npt -pp npt -po npt
gmx mdrun -v -deffnm npt

gmx grompp -f input/pro.mdp -c npt.gro -p topol.top -o pro -pp pro -po pro
gmx mdrun -v -deffnm pro
