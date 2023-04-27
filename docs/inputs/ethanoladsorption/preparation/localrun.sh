#!/bin/sh

gmx grompp -f inputs/min.mdp -c solvated_vacuum.gro -p topol.top -o min -pp min -po min -maxwarn 1
gmx mdrun -v -deffnm min

gmx grompp -f inputs/nvt.mdp -c min.gro -p topol.top -o nvt -pp nvt -po nvt -maxwarn 1
gmx mdrun -v -deffnm nvt


#/work/sgravelle/Softwares/gromacs-install/bin/gmx grompp -f inputs/pro.mdp -c npt.gro -p topol.top -o pro -pp pro -po pro -maxwarn 1
#/work/sgravelle/Softwares/gromacs-install/bin/gmx mdrun -v -deffnm pro
