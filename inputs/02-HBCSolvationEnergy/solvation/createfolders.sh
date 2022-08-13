#/bin/bash
# delete runall.sh if it exist, then re-create it
if test -f "runall.sh"; then
	rm runall.sh
fi
touch runall.sh
echo '#/bin/bash' >> runall.sh
echo '' >> runall.sh
# folder for analysis
mkdir -p dhdl
# loop on the 21 lambda state
for state in $(seq 0 20); 
do
    # create folder
    DIRNAME=lambdastate${state}
    mkdir -p $DIRNAME
    # copy the topology, inputs, and configuration file in the folder
    cp -r topol.top $DIRNAME
    cp -r ../preparation/npt.gro $DIRNAME/preparedstate.gro
    cp -r inputs $DIRNAME
    # replace the lambda state in both npt_bis and production mdp file
    newline='init-lambda-state = '$state
    linetoreplace=$(cat $DIRNAME/inputs/npt_bis.mdp | grep init-lambda-state)
    sed -i '/'"$linetoreplace"'/c\'"$newline" $DIRNAME/inputs/npt_bis.mdp
    sed -i '/'"$linetoreplace"'/c\'"$newline" $DIRNAME/inputs/pro.mdp
    # create a bash file to launch all the simulations
    echo 'cd '$DIRNAME >> runall.sh
    echo '/work/sgravelle/Softwares/gromacs-install/bin/gmx grompp -f inputs/npt_bis.mdp -c preparedstate.gro -p topol.top -o npt_bis -pp npt_bis -po npt_bis -maxwarn 1' >> runall.sh
    echo '/work/sgravelle/Softwares/gromacs-install/bin/gmx mdrun -v -deffnm npt_bis' >> runall.sh
    echo '/work/sgravelle/Softwares/gromacs-install/bin/gmx grompp -f inputs/pro.mdp -c npt_bis.gro -p topol.top -o pro -pp pro -po pro -maxwarn 1' >> runall.sh
    echo '/work/sgravelle/Softwares/gromacs-install/bin/gmx mdrun -v -deffnm pro' >> runall.sh
    echo 'cd ..' >> runall.sh
    echo '' >> runall.sh
    # create links for the analysis
    cd dhdl
    ln -sf ../$DIRNAME/pro.xvg md$state.xvg
    cd ..    
done

