.. _demo_RigiFlex:

Demo RigiFlex
==========================

Run it by typing ``MMMx demo_RigiFlex``. This example generates logfile ``demo_RigiFlex.log``. The logfile is automatically displayed after completion by the ``#report`` directive at the end.
Note that this demo runs for a long time and creates a sizeable ensemble. Therefore, it is contained in a separate subdirectory of ``MMMx\example_set``.
You may want to run and understand ``demo_Rigi`` and ``demo_FlexRNA`` before and run the ``demo_RigiFlex`` over night. 

Modules used
---------------------------------

Modules ``Prepare``, ``Rigi``, ``FlexRNA``, and ``YasaraRefine`` are used.

Functionality
---------------------------------

This demo is a pipeline that results from concatenating ``demo_Rigi`` and ``demo_FlexRNA``. Please consult the documentation of these two parts for a description of functionality.

The only change required for concatenation concerns the first ``FlexRNA`` call, where the input statement ``addpdb RsmE_RsmZ_rigid_bodies`` is replaced by ``expand RsmE_RsmZ_rba``.
The new input statement expands the rigid-body output of ``Rigi`` to PDB files. The ``FlexRNA`` calls are written in a way that they operate on the whole ensemble.
Ensemble size reduces between these calls, since some of the rigid-body arrangements cannot be connected by (all) RNA linkers.

Outputs
---------------------------------

Outputs are the rigid-body template file ``RsmE_RsmZ_rigid_bodies.pdb``, PDB files of all intermediate models, and the rigid-body arrangement file ``RsmE_RsmZ_rba.mat``.
