.. _source_installation:

Source-code version
==========================

Requirements
---------------------------------

The source-code version of MMMx runs under Matlab 2019b or later and requires the Matlab Optimization, 
Global Optimization, Parallel Computing, and Statistics and Machine Learning toolboxes. Please use
the :ref:`executable version<executable_installation>` if you do not have a license for this software, but note that the executbale version is only infrequently updated.

For certain functionality, MMMx relies on :ref:`third-party software<third_party>`, such as `SCWRL4 <http://dunbrack.fccc.edu/SCWRL3.php/>`_,
`ATSAS <https://www.embl-hamburg.de/biosaxs/software.html>`_, `MUSCLE <http://www.drive5.com/muscle/downloads.htm>`_,
`MSMS <http://mgl.scripps.edu/people/sanner/html/msms_home.html>`_, `DSSP <https://swift.cmbi.umcn.nl/gv/dssp/HTML/distrib.html>`_,
`HYDROPRO <http://leonardo.inf.um.es/macromol/programs/hydropro/hydropro.htm>`_, or `YASARA <http://www.yasara.org/>`_.
Except for YASARA structure, which we use for structure refinement, this software is freely available for academic purposes,
but may require registration and a license agreement. 

Download
---------------------------------

The current version of the source code is available at `GitHub <https://github.com/gjeschke/MMMx>`_. Accessing it through GitHub Desktop is the most convenient way of keeping it updated. Otherwise, the package can also be downloaded via the <> Code button and the Downlaod ZIP menu item. Unpack the ZIP file.

Matlab path
---------------------------------

The main MMMx directory and all its subdirectories must be on the Matlab path. Likewise, the directories of all 
required third-party software, often including subdirectories, must be on the Matlab path.

Getting started
---------------------------------

Type ``MMMx help`` or ``MMMx docs``  to access the documentation.

Type ``MMMx demo_distributions.mcx`` to run a first example. Examples exist for all modules.

Type just ``MMMx`` to open a file browser. You can select any MMMx control file (.mcx) to run a modeling pipeline.

For programmatic access from within Matlab, refer to the programmer's documentation.