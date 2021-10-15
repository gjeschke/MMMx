.. _prepare:

Prepare
==========================

This module helps in preparing input PDB files for other MMMx modules. It can cut and merge existing structures, 
reorient structures based on symmetry or on constructing a coarse-grained lipid bilayer model (requires `MSMS <http://mgl.scripps.edu/people/sanner/html/msms_home.html>`_),
add, repair or modify sidechains (requires `SCWRL4 <http://dunbrack.fccc.edu/SCWRL3.php/>`_), superimpose structures (may require `MUSCLE <http://www.drive5.com/muscle/downloads.htm>`_ if sequence-based),
amd renumber residues.

Features are demonstrated in examples ``demo_Prepare.mcx`` and ``demo_RBreference.mcx``. 
The latter example is a two-module pipeline that uses module ``Prepare`` before module :ref:`ExperimentDesign <experiment_design>`. 

The following keywords are supported: