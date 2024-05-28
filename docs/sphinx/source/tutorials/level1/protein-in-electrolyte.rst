.. protein-electrolyte-label:

Protein in electrolyte
**********************

.. container:: hatnote

    Simulating a solvated protein

.. container:: justify

    The goal of this tutorial is to use GROMACS and perform a simple
    molecular dynamics simulation of a protein solvated in an electrolyte. The
    protein will be downloaded from the Protein Data Bank (PDB) :cite:`bank1971protein`
    and solvated in an electrolyte.

.. include:: ../../non-tutorials/needhelp.rst
.. include:: ../../non-tutorials/recommand-salt.rst

Convert the pdb file
====================

.. container:: justify

    Download the *.pdb* file from the |ProteinDataBank|,
    or simply click |1cta.pdb|. The protein is a calcium-binding peptide from site III
    of chicken troponin-C that has been determined using 1H-NMR
    spectroscopy :cite:`shaw1992determination`.

.. |ProteinDataBank| raw:: html

    <a href="https://www.rcsb.org/structure/1CTA" target="_blank">Protein Data Bank</a>

.. |1cta.pdb| raw:: html

    <a href="https://raw.githubusercontent.com/gromacstutorials/gromacstutorials-inputs/main/level1/protein-in-electrolyte/1cta.pdb" target="_blank">here</a>

.. container:: justify

    We first need to create the *.gro* file, i.e. a GROMACS structure file,
    from the *.pdb* file. This can be done using *gmx trjconv*:

.. container:: bash

    gmx trjconv -f 1cta.pdb -s 1cta.pdb -o 1cta.gro -center -box 5 5 5

.. container:: justify

    Choose the group *System* for the centering, and the group *System* as well 
    for the output. A file named *1cta.gro* is created. The generated *.gro*
    file contains 666 atoms, each atom corresponding to one line:

.. container:: bash

    TROPONIN C SITE III - SITE III HOMODIMER
    666
        0ACE      C    1   2.662   4.131   2.701
        0ACE      O    2   2.714   4.036   2.646
        0ACE    CH3    3   2.651   4.147   2.853
    (...)
   35NH2    HN2  664   2.417   3.671   3.192
   69CA      CA  665   3.016   2.279   1.785
   70CA      CA  666   1.859   2.046   1.838
   5.00000   5.00000   5.00000

.. container:: justify

    The last line is the box dimensions in nanometer, which was requested 
    in the *gmx trjconv* command by the *-box 5 5 5* option. All the options
    of *trjconv* can be found in the corresponding page of the GROMACS
    |trjconv-documentation|.

.. |trjconv-documentation| raw:: html

    <a href="https://manual.gromacs.org/current/onlinehelp/gmx-trjconv.html" target="_blank">documentation</a>

.. include:: ../../non-tutorials/accessfile.rst
