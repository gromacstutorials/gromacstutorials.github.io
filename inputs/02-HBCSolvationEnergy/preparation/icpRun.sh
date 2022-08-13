#!/bin/sh

/work/sgravelle/Softwares/gromacs-install/bin/gmx grompp -f inputs/min.mdp -c solvated.gro -p topol.top -o min -pp min -po min -maxwarn 1
/work/sgravelle/Softwares/gromacs-install/bin/gmx mdrun -v -deffnm min

/work/sgravelle/Softwares/gromacs-install/bin/gmx grompp -f inputs/nvt.mdp -c min.gro -p topol.top -o nvt -pp nvt -po nvt -maxwarn 1
/work/sgravelle/Softwares/gromacs-install/bin/gmx mdrun -v -deffnm nvt

/work/sgravelle/Softwares/gromacs-install/bin/gmx grompp -f inputs/npt.mdp -c nvt.gro -p topol.top -o npt -pp npt -po npt -maxwarn 1
/work/sgravelle/Softwares/gromacs-install/bin/gmx mdrun -v -deffnm npt

