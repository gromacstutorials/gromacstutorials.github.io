#!/bin/sh

gmx grompp -f ../inputs/min.mdp -c conf.gro -p topol.top -o min -pp min -po min
gmx mdrun -v -deffnm min

gmx grompp -f ../inputs/nvt.mdp -c min.gro -p topol.top -o nvt -pp nvt -po nvt -maxwarn 1
gmx mdrun -v -deffnm nvt

gmx grompp -f ../inputs/npt.mdp -c nvt.gro -p topol.top -o npt -pp npt -po npt -maxwarn 2
gmx mdrun -v -deffnm npt

gmx grompp -f inputs/nvt.mdp -c npt.gro -p topol.top -o nvt -pp nvt -po nvt -maxwarn 2
gmx mdrun -v -deffnm nvt

gmx grompp -f inputs/video.mdp -c nvt.gro -p topol.top -o video -pp video -po video -maxwarn 2
gmx mdrun -v -deffnm video

gmx trjconv -s video.tpr -f video.xtc -o video.xtc -pbc nojump
