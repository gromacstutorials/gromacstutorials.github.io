.. _contact-before-you-start-label:

Before you start
****************

..  container:: justify

    *GROMACS tutorials* is made of several tutorials that are
    ordered by increasing difficulty.

Required software
=================

..  container:: justify

    The 2024.2 version of GROMACS is required
    to follow the tutorials.

GROMACS (2024.2)
----------------

..  container:: justify

    Download and install the 2024.2 version of GROMACS by following the
    instructions of the |gromacs-manual|.
    Depending on your operative system (i.e. Linux, macOS, or Windows),
    the procedure may differ.

.. |gromacs-manual| raw:: html

   <a href="https://manual.gromacs.org/current/index.html" target="_blank">GROMACS manual</a>

VMD (optional)
--------------

..  container:: justify

    In order to visualize the simulation, the version
    1.9.3 of |VMD| will be used :cite:`humphreyVMDVisualMolecular1996`.
    Some basic instructions for VMD are given on *lammpstutorials*, see
    |VMD-lammps-tutorials|. If you prefer, feel free to use an alternative visualization
    software like |Ovito|.

.. |VMD-lammps-tutorials| raw:: html

   <a href="https://lammpstutorials.github.io/sphinx/build/html/tutorials/vmd/vmd-tutorial.html" target="_blank">this link</a>

.. |VMD| raw:: html

   <a href="https://www.ks.uiuc.edu/Research/vmd" target="_blank">VMD</a>
    
.. |Ovito| raw:: html

   <a href="https://www.ovito.org" target="_blank">Ovito</a>
    
Python (optional)
-----------------

..  container:: justify

    The version 2.6.1 of MDAnalysis is used
    together with the version 3.11.4
    of Python :cite:`vanrossumPythonTutorial, michaud-agrawalMDAnalysisToolkitAnalysis2011`.

..  container:: justify

    To plot the results from the simulations,
    the version 3.5.2 of |Matplotlib Pyplot| is used.

.. |Matplotlib Pyplot| raw:: html

   <a href="https://matplotlib.org/3.5.3/api/_as_gen/matplotlib.pyplot.html" target="_blank">Matplotlib Pyplot</a>

Text editing software
---------------------

..  container:: justify

    To write and edit GROMACS input files, a text editor is required.
    Any text editor will do, such as |gedit|,
    |vim|,
    or |vscode|.
    
.. |gedit| raw:: html

   <a href="https://help.gnome.org/users/gedit/stable/" target="_blank">gedit</a>
    
.. |vim| raw:: html

   <a href="https://www.vim.org/" target="_blank">vim</a>
    
.. |vscode| raw:: html

   <a href="https://code.visualstudio.com/" target="_blank">vscode</a>
    
Find the input scripts
======================

.. include:: ../non-tutorials/accessfile.rst

Recommended reading
===================

..  container:: justify

   To better understand molecular dynamics simulations, I recommend the reading
   of *Understanding molecular simulation* by Daan Frenkel and Berend
   Smit :cite:`frenkelUnderstandingMolecularSimulation2002`, as well as
   *Computer simulation of liquids* by Michael Allen and Dominic Tildesley
   :cite:`allenComputerSimulationLiquids2017`. To understand the basic concepts 
   of fluid and Soft Matter systems, I recommend reading *Basic concepts for
   simple and complex liquids* by Jean-Louis Barrat and Jean-Pierre Hansen
   :cite:`barratBasicConceptsSimple2003`,
   as well as *Theory of simple liquids: with applications to soft matter*
   by Jean-Pierre Hansen and Ian Ranald McDonald :cite:`hansenTheorySimpleLiquids2013a`.