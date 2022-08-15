#!/bin/sh

gmx grompp -f inputs/pull.mdp -c pull.gro -p topol.top -o pull -pp pull -po pull -maxwarn 2 -n index.ndx 
gmx mdrun -v -deffnm pull




