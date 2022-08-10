#!/bin/sh

gmx grompp -f input/min.mdp -c conf.gro -p topol.top -o min -pp min
gmx mdrun -v -deffnm min

gmx grompp -f input/nvt.mdp -c min.gro -p topol.top -o nvt -pp nvt
gmx mdrun -v -deffnm nvt

gmx grompp -f input/npt.mdp -c nvt.gro -p topol.top -o npt -pp npt
gmx mdrun -v -deffnm npt

gmx grompp -f input/run.mdp -c npt.gro -p topol.top -o run -pp run
gmx mdrun -v -deffnm run
