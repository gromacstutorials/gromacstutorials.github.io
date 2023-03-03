.. _bulk-solution-label:

Bulk salt solution
******************

.. container:: hatnote

    The very basics of GROMACS through a
    simple example: a bulk solution of
    SO\ :sub:`4`\ :sup:`2-` and Na\ :sup:`+`.

.. figure:: figures/bulksolution/binary_LJ_fluid.webp
    :alt: Water solution of SO\ :sub:`4`\ :sup:`2-` and Na\ :sup:`+` ions visualized with VMD
    :height: 250
    :align: right

..  container:: justify

    **The objective of this tutorial** is to use the
    open-source code GROMACS to perform a simple molecular
    dynamics simulation: a liquid solution of water mixed
    with sodium (Na\ :sup:`+`) and sulfate
    (SO\ :sub:`4`\ :sup:`2-`) ions. This tutorial
    illustrates several major ingredients of molecular
    dynamics simulations, such as energy minimization,
    thermostating, NVT and NPT equilibration, and
    trajectory visualisation.

    There are **no prerequisite** to follow this tutorial,
    but your life will be easier if you are familiar with
    using the terminal.

.. include:: ../../contact/needhelp.rst

Required softwares
==================

..  container:: justify

    GROMACS must be installed on your machine. You can
    install it following the instructions of the `GROMACS
    website <https://manual.gromacs.org/current/index.html>`__.
    Alternatively, if you are using Ubuntu OS, you can
    simply execute the following command in a terminal:

..  code-block:: bw

    sudo apt-get install gromacs
   
..  container:: justify

    You can verify that GROMACS is indeed installed on your
    computer by typing in a terminal :

..  code-block:: bw

    gmx

..  container:: justify

    You should see the version of GROMACS that has been
    installed. On my computer I see

..  code-block:: bw

    :-) GROMACS - gmx, 2023 (-:

..  container:: justify

as well as a quote (at the bottom), such as

..  code-block:: bw

    GROMACS reminds you: "Computers are like humans - they do everything except think." (John von Neumann)

..  container:: justify

    In addition to GROMACS, you will also need |(1) a basic text editing software|
    such as Vim, Gedit, or Notepad++, |(2) a visualization software|, here I
    will use VMD (note: VMD is free but you have to register to
    the uiuc website in order to download it. If you don't want
    to, you can also use Ovito.), |(3) a plotting tool| like
    XmGrace or pyplot.

.. |LAMMPS website| raw:: html

   <a href="https://lammps.sandia.gov" target="_blank">LAMMPS website</a>

.. |(1) a basic text editing software| raw:: html

   <a href="https://help.gnome.org/users/gedit/stable/" target="_blank">(1) a basic text editing software</a>

.. |(2) a visualization software| raw:: html

   <a href="https://www.ks.uiuc.edu/Research/vmd/" target="_blank">(2) a visualization software</a>

.. |(3) a plotting tool| raw:: html

   <a href="https://plasma-gate.weizmann.ac.il/Grace/" target="_blank">(3) a plotting tool</a>

The input files
===============

..  container:: justify

    In order to run the present simulation using GROMACS,
    we need the 3 following files:

    - 1) a **configuration file** (.gro) containing the
         initial positions of the atoms and the box
         dimensions,
    - 2) a **topology file** (.top) containing
         information about the force field (e.g. interaction
         parameters, molecular constraints),
    - 3) and an **input file** (.mdp) containing the
         parameters of the simulation (e.g. temperature,
         timestep).


The configuration file (.gro)
-----------------------------

..  container:: justify

    For the present simulation, the initial atoms
    positions and box size are given in a conf.gro file
    (Gromos87 format) that you can download by clicking
    `here <https://raw.githubusercontent.com/gromacstutorials/gromacstutorials.github.io/main/inputs/01-Na2SO4solution/conf.gro>`__.
    Save the conf.gro file in a folder. The file looks
    like that:

..  code-block:: bw

    Na2SO4 solution
    3557
        1  SO4   O1    1   2.725   1.999   2.077
        1  SO4   O2    2   2.679   2.091   1.839
        1  SO4   O3    3   2.505   2.127   2.028
        1  SO4   O4    4   2.542   1.890   1.930
        1  SO4    S    5   2.612   2.027   1.968
    (...)
    898  SOL  HW1 3555   3.004   3.021   2.945
    898  SOL  HW2 3556   3.004   2.869   2.945
    898  SOL   MW 3557   2.955   2.945   2.945
    3.10000   3.10000   3.10000

..  container:: justify

    **Description:** The first line is a comment, the
    second is the total number of atoms, and the last
    line is the box dimension in nanometer (nm). Between
    the second and the last line there is one line per
    atom. Each line indicates, from left to right, the
    residue Id (the atoms of the same
    SO\ :sub:`4`\ :sup:`2-` ion have the same residue
    Id), the residue name, the atom name, the atom Id,
    and finally the atom position (x, y, and z
    coordinate in nm).
    
    **Remark:** The format of conf.gro file is fixed,
    all columns are in a fixed position. For example,
    the first five columns are for the residue number.
    This conf.gro file can be visualized using VMD by
    typing in the terminal:

..  code-block:: bash

     vmd conf.gro

..  container:: justify

    You have to play with atoms' representation and color
    to make it look better than it is by default.

#..  code-block:: gromacs
#   :caption: *to be copied in input1.lammps*
