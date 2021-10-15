.. _ensemble_list:

.. highlight:: rst

Ensemble lists
====================

MMMx represents ensemble structures by sets of population-weighted conformers. The standard format is an ensemble list
with extension ``.ens``. 

The file starts with a comment line (initiated by a ``%`` character) the reports the name and the date and time of generation.
This is followed by `C` lines (`C` is the number of conformers in the ensemble) that each contain a PDB file name and a population. 
The populations add up to unity.