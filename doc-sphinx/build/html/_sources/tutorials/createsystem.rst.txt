.. _create-system-label:

Create topology
***************

.. container:: hatnote

    Writing the topology file. 

.. figure:: figures/bulksolution/first-light.png
    :alt: Water solution of SO\ :sub:`4`\ :sup:`2-` and Na\ :sup:`+` ions visualized with VMD
    :class: only-light
    :height: 250
    :align: right

.. figure:: figures/bulksolution/first-dark.png
    :alt: Water solution of SO\ :sub:`4`\ :sup:`2-` and Na\ :sup:`+` ions visualized with VMD
    :class: only-dark
    :height: 250
    :align: right

..  container:: justify

    The objective of this tutorial is to write a 
    simple topology file using python, by placing
    molecules and ions in an empty box. 

    The topology file will be used in ref:`bulk-solution-label`. 
    If you are only interested in learning GROMACS, jump directly
    in ref:`bulk-solution-label`.

.. include:: ../contact/needhelp.rst

What is a .gro file?
====================

A .gro file contains the initial positions of all the atoms 
of a simulation, and can be read by GROMACS. Its structure 
is the following:

..  code-block:: bw

    Name of the system
    number-of-atoms
    residue-number residue-name atom-name atom-number atom-positions (x3) # first atom
    residue-number residue-name atom-name atom-number atom-positions (x3) # second atom
    residue-number residue-name atom-name atom-number atom-positions (x3) # third atom
    (...)
    residue-number residue-name atom-name atom-number atom-positions (x3) # penultimate atom
    residue-number residue-name atom-name atom-number atom-positions (x3) # last atom
    box-size (x3)

..  container:: justify

    One particularity of .gro file format, each column must be located at a fixed position, see |conf.gro-manual|.

.. |conf.gro-manual| raw:: html

    <a href="https://manual.gromacs.org/documentation/current/reference-manual/file-formats.html#gro" target="_blank">the GROMACS manual</a>

Molecule/ions definitions
=========================

Open a blank python script, call it *molecules.py*, and copy the following lines in it:

..  code-block:: python
    :caption: *to be copied in molecules.py*

    import numpy as np

    # define SO4 ion
    def SO4_ion():
        Position = np.array([[0.1238,  0.0587,   0.1119], \
            [0.0778,   0.1501,  -0.1263], \
            [-0.0962,   0.1866,   0.0623], \
            [-0.0592,  -0.0506,  -0.0358],\
            [0.0115,   0.0862,   0.0030]])
        Type = ['OS', 'OS', 'OS', 'OS', 'SO']
        Name = ['O1', 'O2', 'O3', 'O4', 'S1']
        Resname = 'SO4'
        return Position, Type, Resname, Name

    # define Na ion
    def Na_ion():
        Position = np.array([[0, 0, 0]])
        Type = ['Na']
        Name = ['Na1']
        Resname = 'Na'
        return Position, Type, Resname, Name

    # define water molecule
    def H20_molecule():
        Position = np.array([[ 0.    ,  0.    ,  0.    ], \
            [ 0.05858,  0.0757 ,  0.    ], \
            [ 0.05858, -0.0757 ,  0.    ], \
            [ 0.0104 ,  0.    ,  0.    ]])
        Type = ['OW', 'HW', 'HW', 'MW']
        Name = ['OW1', 'HW1', 'HW2', 'MW1']
        Resname = 'Sol'
        return Position, Type, Resname, Name

..  container:: justify

    Each function corresponds to a residue, and contains 
    the positions, type, and name of all the atoms, as well as the name of the residue. 
    These function will be called every time we will need to place a residue in
    our system.

    The sulfide (:math:`\text{SO}_4^{2-}`), sodium (:math:`\text{Na}^{+}`) and water
    molecules (:math:`\text{H}_2\text{O}`) look like that, respectively:

.. figure:: figures/createsystem/molecule-light.png
    :alt: Gromacs tutorial : Initial water molecule, sodium, and sulfide ions.
    :class: only-light

.. figure:: figures/createsystem/molecule-dark.png
    :alt: Gromacs tutorial : Initial water molecule, sodium, and sulfide ions.
    :class: only-dark

    Oxygen atoms are in red, hydrogen atoms in white, sodium atom in blue, and 
    sulfur atom in yellow. The fourth point (MW) of the water molecule is not 
    visible.

Preparing the system
====================

..  container:: justify

    We first need to define the basic parameters, such as 
    the number of residue we want, or the box size, and
    initialize some lists and counters. 

    Next to *molecule.py*, create a nez Python file called
    *generategro.py*, and copy the following lines into it:

..  code-block:: python
    :caption: *to be copied in generategro.py*

    import numpy as np
    from molecules import SO4_ion, Na_ion, H20_molecule

    # define the box size
    Lx, Ly, Lz = [3.6]*3
    box = np.array([Lx, Ly, Lz])

..  container:: justify

    Here box is an array containing the box size along all 3 coordinates of space,
    respectively Lx, Ly, and Lz. Here a cubic box of lateral dimension 3.6 nm is used.

.. include:: ../contact/contactme.rst
