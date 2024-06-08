.. _solvation-energy-label:

Molecule solvation energy
*************************

.. container:: hatnote

    Calculating the free energy of solvation of a graphene-like molecule

.. figure:: ../figures/level3/solvation-energy/video-HBC-light-2.webp
    :alt: HBC (graphene-like) molecule in water
    :class: only-light
    :height: 250
    :align: right

.. figure:: ../figures/level3/solvation-energy/video-HBC-dark-2.webp
    :alt: HBC (graphene-like) molecule in water
    :class: only-dark
    :height: 250
    :align: right

..  container:: justify

    The objective of this tutorial is to use GROMACS
    to perform a molecular simulation of a large molecule in water. 
    By progressively switching off the interactions between the molecule and
    water, the free energy of solvation will be calculated. 

..  container:: justify

    The large and flat molecule used here is a graphene-like and
    discoid molecule named hexabenzocoronene (HBC) of
    formula :math:`\text{C}_{42}\text{H}_{18}`. The TIP4P/epsilon
    model is used for the water :cite:`fuentes2014non`.

.. include:: ../../non-tutorials/recommand-salt.rst
.. include:: ../../non-tutorials/needhelp.rst
.. include:: ../../non-tutorials/GROMACS2024.2.rst

Input files
===========

..  container:: justify

    Create two folders side-by side. Name the two folders *preparation/* and
    *solvation/*, and go to *preparation/*.

..  container:: justify

    Download the configuration files for the HBC molecule
    by clicking |FJEW_allatom_optimised_geometry.pdb|, and save it in
    the *preparation/* folder. The molecule 
    was downloaded from the |atb-HBC| on the automated topology builder
    (ATB) :cite:`malde2011automated`.

.. |FJEW_allatom_optimised_geometry.pdb| raw:: html

    <a href="https://raw.githubusercontent.com/gromacstutorials/gromacstutorials-inputs/main/level3/solvation-energy/preparation/FJEW_allatom_optimised_geometry.pdb" target="_blank">this page</a>

.. |atb-HBC| raw:: html

   <a href="https://atb.uq.edu.au/molecule.py?molid=151371#panel-md" target="_blank">atb</a>

Create the configuration file
-----------------------------

..  container:: justify

    First, let us convert the *.pdb* file into a *.gro* file
    within a cubic box of lateral size 3.5 nanometers using the *gmx trjconv*
    command. Type the following command in a terminal:

..  code-block:: bw

    gmx trjconv -f FJEW_allatom_optimised_geometry.pdb -s FJEW_allatom_optimised_geometry.pdb -o hbc.gro -box 3.5 3.5 3.5 -center  

..  container:: justify

    Select *system* for both centering and output.

.. figure:: ../figures/level3/solvation-energy/hbc-light.png
    :alt: Gromacs initial configuration of HBC graphene molecule
    :class: only-light
    :height: 250
    :align: center

.. figure:: ../figures/level3/solvation-energy/hbc-dark.png
    :alt: Gromacs initial configuration of HBC graphene molecule
    :class: only-dark
    :height: 250
    :align: center

..  container:: figurelegend

    Figure: HBC molecule as seen with VMD with carbon atoms in gray and hydrogen
    atoms in white. The honeycomb structure of the HBC is similar  to that
    of graphene.

..  container:: justify

    Alternatively, you can download the |solvation-hbc.gro| I have generated,
    and continue with the tutorial.

.. |solvation-hbc.gro| raw:: html

    <a href="https://raw.githubusercontent.com/gromacstutorials/gromacstutorials-inputs/main/level3/solvation-energy/preparation/hbc.gro" target="_blank">hbc.gro</a>

Create the topology file
------------------------

..  container:: justify

    Create a folder named *ff/* within the *preparation/* folder.
    Copy the force field parameters from the following *zip* file by 
    clicking |ff-itp.zip|. Both *FJEW_GROMACS_G54A7FF_allatom.itp* file
    and *gromos54a7_atb.ff/* folder were downloaded from the |atb-HBC|.

.. |ff-itp.zip| raw:: html

    <a href="https://raw.githubusercontent.com/gromacstutorials/gromacstutorials-inputs/main/level3/solvation-energy/preparation/ff-itp.zip" target="_blank">here</a>

..  container:: justify

    Then, let us write the topology (*top*) file by simply
    creating a blank file named *topol.top* within the
    *preparation/* folder. Copy the following lines into *topol.top*:

..  code-block:: bw

    #include "ff/gromos54a7_atb.ff/forcefield.itp"
    #include "ff/FJEW_GROMACS_G54A7FF_allatom.itp"

    [ system ]
    Single HBC molecule

    [ molecules ]
    FJEW 1

Solvate the HBC in water
------------------------

..  container:: justify

    Let us add water molecules to the system. First, download the tip4p
    water configuration (*.gro*) file |tip4p.gro|,
    and copy it in the *preparation/* folder. This is a fourth point water
    model with additional massless site where the charge of the
    oxygen atom is placed. Then, in
    order to add tip4p water molecules to both *.gro* and
    *.top* file, use the *gmx solvate* command as follows:

.. |tip4p.gro| raw:: html

    <a href="https://raw.githubusercontent.com/gromacstutorials/gromacstutorials-inputs/main/level3/solvation-energy/preparation/tip4p.gro" target="_blank">here</a>

..  code-block:: bw

    gmx solvate -cs tip4p.gro -cp hbc.gro -o solvated.gro -p topol.top

..  container:: justify

    The new *solvated.gro* file contains all *8804* atoms from the HBC
    molecule (called FJEW) and the water molecules. A new line
    *SOL 2186* also appeared in the topology *.top* file:

..  code-block:: bw

    [ molecules ]
    FJEW 1
    SOL 2186

..  container:: justify

    Alternatively, you can download the *solvated.gro* file I have generated by
    clicking |solvated.gro|, and continue with the tutorial.

.. |solvated.gro| raw:: html

    <a href="https://raw.githubusercontent.com/gromacstutorials/gromacstutorials-inputs/main/level3/solvation-energy/preparation/solvated.gro" target="_blank">here</a>

..  container:: justify

    Finally, save the topology file for the water, the |h2o.itp| file, in
    the *ff/* folder and add the *#include "ff/h2o.itp"* line to the *topol.top*
    file:

.. |h2o.itp| raw:: html

    <a href="https://raw.githubusercontent.com/gromacstutorials/gromacstutorials-inputs/main/level3/solvation-energy/preparation/ff/h2o.itp" target="_blank">h2o.itp</a>

..  code-block:: bw

    #include "ff/gromos54a7_atb.ff/forcefield.itp"
    #include "ff/FJEW_GROMACS_G54A7FF_allatom.itp"
    #include "ff/h2o.itp"

System equilibration
====================

..  container:: justify

    The system is now ready for the simulations. Let us first equilibrate it
    before measuring the solvation energy of the HBC molecule. 

Energy minimization
-------------------

..  container:: justify

    Create an *inputs/* folder inside the *preparation/* folder,
    and create a new blank file called
    *min.mdp* in it. Copy the following lines into *min.mdp*:

..  code-block:: bw

    integrator = steep
    nsteps = 50

    nstxout = 10

    cutoff-scheme = Verlet
    nstlist = 10
    ns_type = grid

    couple-intramol=yes

    vdw-type = Cut-off
    rvdw = 1.2

    coulombtype = pme
    fourierspacing = 0.1
    pme-order = 4
    rcoulomb = 1.2

    constraint-algorithm = lincs
    constraints = hbonds

..  container:: justify

    All these lines have been seen in the previous
    tutorials. In short, with this input script, GROMACS will perform a
    steepest descent by updating the atom positions
    according to the largest forces directions, until
    the energy and maximum force reach a reasonable value. 
    
..  container:: justify

    Apply the minimization to the solvated box using *gmx grompp* and *gmx mdrun*:

..  code-block:: bash

    gmx grompp -f inputs/min.mdp -c solvated.gro -p topol.top -o min -pp min -po min -maxwarn 1
    gmx mdrun -v -deffnm min

..  container:: justify

    Here, the *-maxwarn 1* option is used to ignore a WARNING from GROMACS
    about some issue with the force field. For this tutorial, we can safely
    ignore this WARNING.

..  container:: justify

    Let us visualize the atom trajectories during the
    minimization step using VMD by typing:

..  code-block:: bash

    vmd solvate.gro min.trr

.. figure:: ../figures/level3/solvation-energy/minimize-light.png
    :alt: Gromacs HBC (graphene) molecule after minimization in water
    :class: only-light
    :height: 400
    :align: center

.. figure:: ../figures/level3/solvation-energy/minimize-dark.png
    :alt: Gromacs HBC (graphene) molecule after minimization in water
    :class: only-dark
    :height: 400
    :align: center

.. container:: figurelegend

    Figure: The system after energy minimization showing the HBC molecule in water.
    The water is represented as a transparent field. 

NVT and NPT equilibration 
-------------------------

..  container:: justify

    Let us perform successively a NVT and a NPT relaxation steps. Copy the |solvation-nvt.mdp|
    and the |solvation-npt.mdp| files into the inputs folder, and run them both using:

.. |solvation-nvt.mdp| raw:: html

    <a href="https://raw.githubusercontent.com/gromacstutorials/gromacstutorials-inputs/main/level3/solvation-energy/preparation/inputs/nvt.mdp" target="_blank">nvt.mdp</a>

.. |solvation-npt.mdp| raw:: html

    <a href="https://raw.githubusercontent.com/gromacstutorials/gromacstutorials-inputs/main/level3/solvation-energy/preparation/inputs/npt.mdp" target="_blank">npt.mdp</a>

..  code-block:: bash

    gmx grompp -f inputs/nvt.mdp -c min.gro -p topol.top -o nvt -pp nvt -po nvt -maxwarn 1
    gmx mdrun -v -deffnm nvt
    gmx grompp -f inputs/npt.mdp -c nvt.gro -p topol.top -o npt -pp npt -po npt -maxwarn 1
    gmx mdrun -v -deffnm npt

.. figure:: ../figures/level3/solvation-energy/nvtnpt-light.webp
    :alt: Gromacs HBC (graphene) molecule during NVT and NPT equilibration
    :class: only-light
    :height: 400
    :align: center

.. figure:: ../figures/level3/solvation-energy/nvtnpt-dark.webp
    :alt: Gromacs HBC (graphene) molecule during NVT and NPT equilibration
    :class: only-dark
    :height: 400
    :align: center

.. container:: figurelegend

    Figure: Movie showing the motion of the atoms during the
    NVT and NPT equilibration steps. For clarity, the
    water molecules are represented as a continuum field.

Solvation energy measurement
============================

..  container:: justify

    The equilibration of the system is complete. Let us perform the solvation
    free energy calculation, for which 21 independent
    simulations will be performed.
     
..  admonition:: About free energy calculation
    :class: info
     
    The interactions between the HBC molecule and the water
    are progressively turned-off, thus effectively
    mimicking the HBC molecule moving from bulk water
    to vacuum. Then, the free energy difference between the fully solvated
    and the fully decoupled configurations is measured.

..  container:: justify

    Within the *solvation/* folder, create an *inputs/*
    folders. Copy the two following
    |solvation-npt-bis.mdp| and |solvation-pro.mdp| files in it.

.. |solvation-npt-bis.mdp| raw:: html

    <a href="https://raw.githubusercontent.com/gromacstutorials/gromacstutorials-inputs/main/level3/solvation-energy/solvation/inputs/npt_bis.mdp" target="_blank">npt_bis.mdp</a>

.. |solvation-pro.mdp| raw:: html

    <a href="https://raw.githubusercontent.com/gromacstutorials/gromacstutorials-inputs/main/level3/solvation-energy/solvation/inputs/pro.mdp" target="_blank">pro.mdp</a>

..  container:: justify

    Both files contain the following commands that are
    related to the free energy calculation:

..  code-block:: bw

    free_energy = yes
    vdw-lambdas = 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
    coul-lambdas = 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0
    sc-alpha = 0.5
    sc-power = 1
    init-lambda-state =  0
    couple-lambda0 = none
    couple-lambda1 = vdw-q
    nstdhdl = 100
    calc_lambda_neighbors = -1
    couple-moltype = FJEW

..  container:: justify

    These lines specify that the decoupling between the
    molecule of interest (here FJEW) and the rest of
    the system (here water) must be done by
    progressively turning off van der Waals and Coulomb
    interactions. The parameter *nstdhdl* controls the frequency at
    which information are being printed in a xvg file during
    the production run.

..  container:: justify

    In addition, the stochastic integrator 'sd' is used
    instead of 'md', as it provides a better sampling,
    which is crucial here, particularly when the HBC
    and the water molecules are not coupled.

..  container:: justify

    Copy as well the following |solvation-topol.top| file within the
    *solvation/* folder (the only difference with the previous one if the path
    to the *ff/* folder).

.. |solvation-topol.top| raw:: html

    <a href="https://raw.githubusercontent.com/gromacstutorials/gromacstutorials-inputs/main/level3/solvation-energy/solvation/topol.top" target="_blank">topol.top</a>

..  container:: justify

    We need to create 21 folders, each containing the
    input files with different value of
    init-lambda-state (from 0 to 21). To do so, create
    a new bash file fine within the 'solvation/'
    folder, call it 'createfolders.sh' can copy the
    following lines in it:

..  code-block:: bw

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
        echo 'gmx grompp -f inputs/npt_bis.mdp -c preparedstate.gro -p topol.top -o npt_bis -pp npt_bis -po npt_bis -maxwarn 1' >> runall.sh
        echo 'gmx mdrun -v -deffnm npt_bis -nt 4' >> runall.sh
        echo 'gmx grompp -f inputs/pro.mdp -c npt_bis.gro -p topol.top -o pro -pp pro -po pro -maxwarn 1' >> runall.sh
        echo 'gmx mdrun -v -deffnm pro -nt 4' >> runall.sh
        echo 'cd ..' >> runall.sh
        echo '' >> runall.sh
        # create links for the analysis
        cd dhdl
        ln -sf ../$DIRNAME/pro.xvg md$state.xvg
        cd ..    
    done

..  container:: justify

    Change the *-nt 4* to use a different number of thread if necessary
    or/and possible.

..  container:: justify

    Execute the bash script by typing:

..  code-block:: bash

    bash createfolders.sh

..  container:: justify

    The bash file creates 21 folders, each containing
    the input files with init-lambda-state from 0 to
    21, as well as a 'topol.top' file and a
    'preparedstate.gro' corresponding to the last
    state of the system simulated in the
    'preparation/' folder. Run all 21 simulations by executing the 'runall.sh' script:

..  code-block:: bash

    bash runall.sh

..  container:: justify

    This may take a while, depending on your computer.

    When the simulation is complete, go the dhdl folder, and type:

..  code-block:: bash

    gmx bar -f *.xvg

..  container:: justify

    The value of the solvation energy is printed in the terminal:

..  code-block:: bash

    total 0 - 20, DG -37.0 +/- 8.40

..  container:: justify

    The present simulations are too short to give a
    reliable result. To accurately measure the solvation
    energy of a molecule, use much longer equilibration
    (typically one nanosecond) and production runs
    (typically several nanoseconds).

.. include:: ../../non-tutorials/accessfile.rst
