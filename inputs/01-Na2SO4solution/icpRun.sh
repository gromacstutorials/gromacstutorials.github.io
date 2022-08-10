#!/bin/sh

/work/sgravelle/Softwares/gromacs-install/bin/gmx grompp -f input/min.mdp -c conf.gro -p topol.top -o min -pp min -po min
/work/sgravelle/Softwares/gromacs-install/bin/gmx mdrun -v -deffnm min

/work/sgravelle/Softwares/gromacs-install/bin/gmx grompp -f input/nvt.mdp -c min.gro -p topol.top -o nvt -pp nvt -po nvt
/work/sgravelle/Softwares/gromacs-install/bin/gmx mdrun -v -deffnm nvt

/work/sgravelle/Softwares/gromacs-install/bin/gmx grompp -f input/npt.mdp -c nvt.gro -p topol.top -o npt -pp npt -po npt
/work/sgravelle/Softwares/gromacs-install/bin/gmx mdrun -v -deffnm npt

/work/sgravelle/Softwares/gromacs-install/bin/gmx grompp -f input/pro.mdp -c npt.gro -p topol.top -o pro -pp pro -po pro
/work/sgravelle/Softwares/gromacs-install/bin/gmx mdrun -v -deffnm pro
