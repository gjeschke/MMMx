.. _ensemble_list:

.. highlight:: rst

Ensemble lists
====================

MMMx represents ensemble structures by sets of weighted conformers. The standard format is an ensemble list
with extension ``.ens``. 

The file starts with a comment line (initiated by a ``%`` character) the reports the name and the date and time of generation.
This is followed by `C` lines (`C` is the number of conformers in the ensemble) that each contain a PDB file name and a weight. 
The weights are considered as populations and thus add up to unity. Ensembles can be archived, by generating a ZIP file that 
contains the ensemble list as well as the individual PDB files of the conformers (see ``archive`` keyword in EnsembleAnalysis. 