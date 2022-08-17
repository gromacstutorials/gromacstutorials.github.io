#!/bin/sh

#/work/sgravelle/Softwares/gromacs-install/bin/gmx grompp -f inputs/min.mdp -c solvated_vacuum.gro -p topol.top -o min -pp min -po min -maxwarn 1
#/work/sgravelle/Softwares/gromacs-install/bin/gmx mdrun -v -deffnm min

#/work/sgravelle/Softwares/gromacs-install/bin/gmx grompp -f inputs/nvt.mdp -c min.gro -p topol.top -o nvt -pp nvt -po nvt -maxwarn 1
#/work/sgravelle/Softwares/gromacs-install/bin/gmx mdrun -v -deffnm nvt

/work/sgravelle/Softwares/gromacs-install/bin/gmx grompp -f inputs/pull.mdp -c ../preparation/nvt.gro -p topol.top -o pull -pp pull -po pull -n index.ndx -maxwarn 1
/work/sgravelle/Softwares/gromacs-install/bin/gmx mdrun -v -deffnm pull

rm *#*

