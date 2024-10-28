.. _executable_installation:

Executable
==========================

Requirements
---------------------------------

The executable version of MMMx runs under Windows. Please contact gjeschke@ethz.ch if you need a version for a different operating system.

MMMx version 1.0 requires the free Matlab Runtime Component 2021a or higher. The installer will offer to add it, if this component is
not yet installed on your system.

For certain functionality, MMMx relies on :ref:`third-party software<third_party>`, such as `SCWRL4 <http://dunbrack.fccc.edu/SCWRL3.php/>`_,
`ATSAS <https://www.embl-hamburg.de/biosaxs/software.html>`_, `MUSCLE <http://www.drive5.com/muscle/downloads.htm>`_,
`MSMS <http://mgl.scripps.edu/people/sanner/html/msms_home.html>`_, `DSSP <https://swift.cmbi.umcn.nl/gv/dssp/HTML/distrib.html>`_,
`HYDROPRO <http://leonardo.inf.um.es/macromol/programs/hydropro/hydropro.htm>`_, or `YASARA <http://www.yasara.org/>`_.
Except for YASARA structure, which we use for structure refinement, this software is freely available for academic purposes,
but may require registration and a license agreement. 

Download
---------------------------------

The Windows installer ``MMMx_Installer_web.exe`` and example library ``MMMx_examples.zip`` need to be downloaded separately. 

Installation
---------------------------------

Double-click the installer and follow the instructions. Unpack the example ZIP file to a directory of your choice.

Windows path
---------------------------------

It is advisable to add the directory, where you installed MMMx, to the Windows path.
Otherwise, you can run MMMx only from this directory. Likewise, the directories of all 
required third-party software must be on the Windows path.

Getting started
---------------------------------

MMMx can be started by double-clicking the executable, 
by double-clicking an icon on the Desktop (if this was added during installation), or from a 
Windows command prompt. In the first two cases, a file browser window opens.
You can select any MMMx control file (.mcx) to run a modeling pipeline.

When running from the command window, the file browser window opens if you just type ``MMMx``.
Type ``MMMx help`` or ``MMMx docs``  to access the documentation.

Run ``MMMx distributions.mcx`` for a first test. If you do this from a Windows command prompt, 
you need to change into the MMMx examples directory before.