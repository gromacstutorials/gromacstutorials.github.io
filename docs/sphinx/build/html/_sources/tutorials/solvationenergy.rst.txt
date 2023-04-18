.. _bulk-solution-label:

Molecule solvation energy
*************************

.. container:: hatnote

    Free energy solvation calculation of a
    disk-like molecule

.. figure:: figures/solvationenergy/banner.png
    :alt: HBC (graphene-like) molecule in water
    :height: 250
    :align: right

..  container:: justify

    **The objective of this tutorial** is to use GROMACS
    to perform a molecular dynamics simulation, and to
    calculate the free energy of solvation of a
    graphene-like molecule named hexabenzocoronene.

Input files
===========

..  container:: justify

    Create two folders named 'preparation/' and
    'solvation' in the same directory. Go to
    'preparation/'.

    Download the configuration files for the HBC molecule
    from the atb repository: click
    `here <https://atb.uq.edu.au/molecule.py?molid=151371#panel-md>`__,
    download the structure file 'All-Atom PDB (optimised
    geometry)' and place it in the 'preparation/' folder.

Create the configuration file
-----------------------------

..  container:: justify

    First, let us convert the pdb file into a gro file
    within a box of finite size using trj conv:

..  code-block:: bw

    gmx trjconv -f FJEW_allatom_optimised_geometry.pdb -s FJEW_allatom_optimised_geometry.pdb -o hbc.gro -box 3 3 3 -center  

..  container:: justify

    Select 'system' for both centering and output. If you
    open the hbc.gro file with VMD, you will see:

.. figure:: figures/solvationenergy/hbc.png
    :alt: GROMACS tutorial : HBC (graphene) molecule with VMD
    :height: 150

    HBC molecule with carbon atoms in gray and hydrogen
    atoms in white.

..  container:: justify

    You can also download the |hbc.gro|
    I have generated.

.. |hbc.gro| raw:: html

    <a href="../../../../inputs/02-HBCSolvationEnergy/preparation/hbc.gro" target="_blank">here</a>

Create the topology file
------------------------

..  container:: justify

    From the `same atb
    page <https://atb.uq.edu.au/molecule.py?molid=151371#panel-md>`__,
    copy the 'GROMACS G54A7FF All-Atom (ITP file)' and
    place it in a folder named 'ff/' and located within
    the 'preparation/' folder. Within 'ff/' download as
    well the GROMACS top file named `Gromacs 4.5.x-5.x.x
    54a7 <https://atb.uq.edu.au/forcefield_files/atb_gromacs/5/gromos54a7_atb.ff.tar.gz>`__
    containing all the force field parameters.
    Then, let us write the topology file by simply
    creating a blank file named 'topol.top' within the
    'preparation/' folder, and copying in it:

..  code-block:: bw

    #include "ff/gromos54a7_atb.ff/forcefield.itp"
    #include "ff/FJEW_GROMACS_G54A7FF_allatom.itp"

    [ system ]
    Single HBC molecule

    [ molecules ]
    FJEW 1

Add the water
-------------

..  container:: justify

    Let us add water molecules. First download the tip4p
    water configuration file |tip4p.gro|,
    and copy it in the 'preparation/' folder. Then, in
    order to add (tip4p) water molecules to both gro and
    top file, use the gmx solvate command as follow:

..  code-block:: bw

    gmx solvate -cs tip4p.gro -cp hbc.gro -o solvated.gro -p topol.top

.. |tip4p.gro| raw:: html

    <a href="../../../../inputs/02-HBCSolvationEnergy/preparation/tip4p.gro" target="_blank">here</a>

..  container:: justify

    You should see the following message:

..  code-block:: bw

    Processing topology
    Adding line for 887 solvent molecules with resname (SOL) to topology file (topol.top)

..  container:: justify

    and a new line 'SOL 887' in the topology file:

..  code-block:: bw

    [ molecules ]
    FJEW 1
    SOL 887

.. include:: ../../contact/recommand-lj.rst

.. include:: ../../contact/needhelp.rst

.. include:: ../../contact/contactme.rst
