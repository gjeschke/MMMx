.. _demo_Rigi:

Demo Rigi
==========================

This example generates logfile ``demo_rigi.log``. The logfile is automatically displayed after completion by the ``#report`` directive at the end.

Modules used
---------------------------------

Modules ``Prepare`` and ``Rigi`` are used.

Source code
------------

.. literalinclude :: ..\..\example_set\demo_Rigi.mcx
   :language: matlabsession


Functionality
---------------------------------

This example reassembles the RsmE/RsmZ complex, excluding flexible RNA linkers, from its parts and from simulated distance distribution restraints.   

First, the prepare module is used to generate a rigid-body template containing only the copies of the RsmE protein and the RNA stemloops of RsmZ.
To this end, structure ``2MF1`` is downloaded from the PDB and stored in entity ``RSM0``. The six RsmE copies are chains A-F, residues 1-59. They are extracted from the first model ``{1}`` of this NMR structure ensemble.
The stemloops are extracted from chain G of this model. All these parts are merged into entity ``RSM1``. This entity is then saved to file ``Rsme_Rsmz_rigid_bodies.pbd`` with pseudo-PDB identifier ``RSM1``.

In the rigi module, this PDB file is specified as rigid-body template. The rigid bodies are separated from each other. A maximum number of ``50000`` trials is specified.
The rigid-body arrangements are saved in MMMx-internal format to file ``RsmE_RsmZ_rba.mat``. The rigid-body arrangements cover 50\% of total probability.

Rigid body 1 is assmebled from an RsmE dimer (chains A and B of the tempalte) and the stemlops attached to this dimer (chains G and H of the template, corresponding to (G)20-35 and (G)44-56 in PDB structure 2MF1).
Residues 8 and 40 in chain A and residue 40 in chain B are specified as reference labelling sites with spin label ``mtsl``.
Rigid body 2 is assembled analogously. Rigid body 3 features only one stemloop.

Then, 18 distance distribution restraints with mtsl labels are specified between reference sites (core restraints) and six further restraints are specified that involve equivalent sites, which are not reference sites (auxiliary restraints).
Each restraint is specified by the mean distance and standard deviation of the distance distribution.

Finally, the nucleotide links between the stemloops are specified, each by the number of nucleotide links and the maximum distance. For the maximum distance, 6 Å per nucleotide link are assumed. 

Outputs
---------------------------------

Outputs are the rigid-body template file ``RsmE_RsmZ_rigid_bodies.pdb`` and the rigid-body arrangement file ``RsmE_RsmZ_rba.mat``.