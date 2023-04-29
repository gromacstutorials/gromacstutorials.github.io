#!/bin/bash

set -e

for ((i = 0 ; i < 30 ; i++)); do
	x=$(echo "0.05*$(($i+1))" | bc);
	sed 's/to_be_replaced/'$x'/g' inputs/min.mdp > min.mdp
	gmx grompp -f min.mdp -c ../preparation/nvt_1ns.gro -p topol.top -o min.$i -pp min.$i -po min.$i -maxwarn 1 -n index.ndx
	gmx mdrun -v -deffnm min.$i

done

#gmx grompp -f inputs/min.mdp -c ../preparation/nvt_1ns.gro -p topol.top -o min -pp min -po min -maxwarn 2 -n index.ndx 
#gmx mdrun -v -deffnm min

#gmx grompp -f inputs/pull.mdp -c min.gro -p topol.top -o pull -pp pull -po pull -maxwarn 2 -n index.ndx 
#gmx mdrun -v -deffnm pull




