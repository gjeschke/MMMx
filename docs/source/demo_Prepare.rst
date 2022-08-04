.. _demo_Prepare:

Demo Prepare
==========================

Run it by typing ``MMMx demo_Prepare``. This example generates logfile ``demo_Prepare.log``.

Modules used
---------------------------------

Only module ``Prepare`` is used.

Source code
------------

.. literalinclude :: ..\..\example_set\demo_Prepare.mcx
   :language: matlabsession


Functionality
---------------------------------

The first example loads the ensemble structure ``2ADC`` from the PDB into entity ``NMR``. Then it saves a local PDB file ``rb34.pdb`` that contains only residues 334-528 of chain A in only model 3 of the ensemble.

The second example generates a chimera from two structures related to vitamin B12 transporter BtuCDF in a frame related to a lipid bilayer and renumbers one of the chains. 

First, the PDB file ``2QI9`` is downloaded as initial model with internal name ``BtuCDF``. Then the coodinates are translated, so that chains A and B (the transmembrane part) are centered.
A pseudo-symmetry axis is computed that relates chains A and B and the coordinates are transformed so that this pseudo-symmetry axis is the z axis.
This model is saved as file ``BtuCDF_centered.pdb``.

Then, selenium amino acids, which were used for phasing the x-ray data, are replaced by their (native) sulfur counterparts. This model is saved as file ``BtuCDF_deselenated.pdb``.

Now, a bilayer geometry is computed in helix ``bundle`` mode, assuming that the z axis (here the pseudo-symmetry axis) is the bilayer normal (mode ``oriented``).
This model is saved as file ``BtuCDF_bilayer_transform.pdb``. Compared to the previous model, only the z coordinates have changed, so that z = 0 corresponds to the bilayer central plane.
The optimal bilayer thickness and the applied coordinate shift are rported in the logfile.

Then, PDB file ``5M29`` is downlaoded as internal entity BtuF_CBI. This structure is the substrate-binding protein BtuF in complex with cobinamid.
The structure contains a cyanide ion and two glycerol molecules that are not required in the chimera. They are removed in the following.
This model is saved as file ``BtuF_CBI_removal.pdb``.

The substrate-binding protein, which is chain A in ``2M29`` is superimposed on the substrate-binding protein in BtuCDF, where it is chain F. Directive ``align`` uses sequence alignment to match residue numbers.
Directive ``backbone`` superimposes backbone atoms. 

Chain F in model BtuCDF is then replaced by chain A in model BtuF_CBI. Because the cofactor cobinamide is a cofactor of chain A of this model, this introduces the substrate molecule into the chimera.
The chimera is saved to file ``BtuCDF_chimera`` with pseudo-PDB code ``BTUX``.

The final part demonstrates renumbering of chain F, with all new residue numbers decreased by 20. This models is saved to file ``BtuCDF_chimera_renumbered.pdb``. 



Outputs
---------------------------------

All outputs, except for the logfile, are PDB files. All intermediate models are stored here. In most application scenarios, it will be sufficient to save the final model.