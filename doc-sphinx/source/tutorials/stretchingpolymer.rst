.. _stretching-polymer-label:

Polymer in water
****************

.. container:: hatnote

    Solvating a small molecule in water before stretching it

.. figure:: figures/stretchingpolymer/main-dark.png
    :alt: peg molecule in water
    :height: 250
    :align: right
    :class: only-dark

.. figure:: figures/stretchingpolymer/main-light.png
    :alt: peg molecule in water
    :height: 250
    :align: right
    :class: only-light

..  container:: justify

    The goal of this tutorial is to use GROMACS and
    create a small hydrophilic polymer (PEG -
    PolyEthylene Glycol) in a reservoir of water. 
    An all-atom description is used, therefore all species considered here
    are made of charged atoms connected by bonds constraints.

    Once the system is created, a constant stretching force will be applied to both
    ends of the polymer, and its length will be measured with time.

    This tutorial was inspired by a very nice |Liese2017| by Liese and coworkers, in which
    they compare MD simulations with force spectroscopy experiments.

.. |Liese2017| raw:: html

    <a href="https://doi.org/10.1021/acsnano.6b07071" target="_blank">publication</a>

.. include:: ../contact/needhelp.rst

PEG molecule in vacuum
======================

..  container:: justify

    Download the *peg.gro* file for the PEG molecule by clicking |download_H2O.data|.

.. |download_H2O.data| raw:: html

   <a href="../../../../inputs/stretchingpolymer/peg.gro" target="_blank">here</a>

..  container:: justify

    Opening *peg.gro* using VMD, one can see that it consists of 
    a rather long polymer chain main of carbon atoms (in gray),
    oxygen atoms (in red), and hydrogen atoms (in white):

.. figure:: figures/stretchingpolymer/light-PEG.png
   :alt: PEG polymer for molecular dynamics simulation in GROMACS
   :class: only-light
   :width: 500

.. figure:: figures/stretchingpolymer/dark-PEG.png
   :alt: PEG polymer for molecular dynamics simulation in GROMACS
   :class: only-dark
   :width: 500

..  container:: justify

    Create a folder named *peg-in-vacuum/*, and copy
    *peg.gro* in it. Next to *peg.gro* create an empty
    file named *topol.top*, and copy the following lines in it:

..  code-block:: bw
    :caption: *to be copied in topol.top*

    [ defaults ]
    ; nbfunc	comb-rule	gen-pairs	fudgeLJ	fudgeQQ
      1         1           no          1.0     1.0

    ; Include forcefield parameters
    #include "ff/charmm35r.itp"
    #include "ff/peg.itp"

    [ system ]
    ; Name
      PEG

    [ molecules ]
    ; Compound        #mols
      PEG             1

..  container:: justify

    Next to *conf.gro* and *topol.top*, create a folder named *ff/*, and copy the
    following files into it: |download_charmm35r.itp| and |download_peg.itp|

.. |download_charmm35r.itp| raw:: html

   <a href="../../../../inputs/stretchingpolymer/ff/charmm35r.itp" target="_blank">charmm35r.itp</a>
   
.. |download_peg.itp| raw:: html

   <a href="../../../../inputs/stretchingpolymer/ff/peg.itp" target="_blank">peg.itp</a>

.. include:: ../contact/contactme.rst
