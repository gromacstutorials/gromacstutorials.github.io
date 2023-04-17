#!/bin/bash

set -e

for ((i = 0 ; i < 30 ; i++)); do
	x=$(echo "0.065*$(($i+1))" | bc);
	sed 's/to_be_replaced/'$x'/g' inputs/min.mdp > min.mdp
	/work/sgravelle/Softwares/gromacs-install/bin/gmx grompp -f min.mdp -c ../preparation/nvt_1ns.gro -p topol.top -o min.$i -pp min.$i -po min.$i -maxwarn 1 -n index.ndx
	/work/sgravelle/Softwares/gromacs-install/bin/gmx mdrun -v -deffnm min.$i
	
	sed 's/to_be_replaced/'$x'/g' inputs/nvt.mdp > nvt.mdp
	/work/sgravelle/Softwares/gromacs-install/bin/gmx grompp -f nvt.mdp -c min.$i.gro -p topol.top -o nvt.$i -pp nvt.$i -po nvt.$i -maxwarn 1 -n index.ndx
	/work/sgravelle/Softwares/gromacs-install/bin/gmx mdrun -v -deffnm nvt.$i
	
	sed 's/to_be_replaced/'$x'/g' inputs/pro.mdp > pro.mdp
	/work/sgravelle/Softwares/gromacs-install/bin/gmx grompp -f pro.mdp -c nvt.$i.gro -p topol.top -o pro.$i -pp pro.$i -po pro.$i -maxwarn 1 -n index.ndx
	/work/sgravelle/Softwares/gromacs-install/bin/gmx mdrun -v -deffnm pro.$i
done


