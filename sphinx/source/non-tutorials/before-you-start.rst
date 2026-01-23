.. _contact-before-you-start-label:
.. include:: links.rst

Before You Start
****************

The GROMACS tutorials in this series are organized by increasing complexity. If
you're completely new to GROMACS, it's strongly recommended to begin with the
first tutorial on a simple :ref:`bulk-solution-label` to get familiar with the
software and workflow.

Required Software
=================

Each tutorial specifies the minimum required version of GROMACS at the top.
While most examples should work with any recent version, using a different
version may lead to unexpected errors or differences in behavior or output.
When possible, try to match the GROMACS version used in the
tutorial when possible.

GROMACS
-------

Download and install GROMACS by following the instructions in the
|gromacs-manual|. The installation procedure may vary depending on your
operating system (e.g., Linux, macOS, or Windows).

VMD
---

In order to visualize the simulation, install Visual Molecular Dynamics (|VMD|)
:cite:`humphreyVMDVisualMolecular1996`. Version 1.9.3 of VMD was used to
create the images in this website, but more recent versions should work as
well. If you prefer, feel free to use an alternative visualization software.
    
Python (optional)
-----------------

The version 2.6.1 of MDAnalysis is used together with version 3.11.4 of
Python :cite:`vanrossumPythonTutorial1995,michaud-agrawalMDAnalysisToolkitAnalysis2011`.
To plot the results from the simulations, version 3.5.2 of |Matplotlib-Pyplot| is used.
All the script used to generate the plots are available from |GitHub-repository|.

Text editing software
---------------------

To write and edit GROMACS input files, a text editor is required.
Any plain text editor will do, such as |gedit|, |vim|, or |vscode|,
|nano|, |sublime|, and |notepadpp|.
    
Find the input scripts
======================

.. include:: ../non-tutorials/accessfile.rst

Recommended reading
===================

To better understand molecular dynamics simulations, I recommend the reading of
*Understanding Molecular Simulation* by Daan Frenkel and Berend Smit
:cite:`frenkelUnderstandingMolecularSimulation2002`, as well as *Computer
Simulation of Liquids* by Michael Allen and Dominic Tildesley
:cite:`allenComputerSimulationLiquids2017`.

To understand the basic concepts of fluid and Soft Matter systems, I recommend
reading *Basic Concepts for Simple and Complex Liquids* by Jean-Louis Barrat and
Jean-Pierre Hansen :cite:`barratBasicConceptsSimple2003`, as well as *Theory of
Simple Liquids: With Applications to Soft Matter* by Jean-Pierre Hansen and Ian
Ranald McDonald :cite:`hansenTheorySimpleLiquids2013a`.
