#/bin/bash

# Create links for dhdl file
mkdir -p dhdl
cd dhdl
    for j in $(seq 0 20); 
    do
        ln -sf ../lambda_$j/grompp_$j.xvg md_$j.xvg
    done
cd ..    