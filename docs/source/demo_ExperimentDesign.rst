.. _demo_ExperimentDesign:

Demo ExperimentDesign
==========================

Run it by typing ``MMMx demo_ExperimentDesign``. This example generates logfile ``demo_ExperimentDesign.log``.

Modules used
---------------------------------

Only module ``ExperimentDesign`` is used.

Source code
------------

.. literalinclude :: ..\..\example_set\demo_ExperimentDesign.mcx
   :language: matlabsession


Functionality
---------------------------------

The demo imports a T4 Lysozyme structure with PDB identifier ``2LZM`` and stores it with internal identifier ``T4L``.
It also imports a RNA recognition motif with bound RNA stemloop with PDB identifier ``2ADC`` and stores it with internal identifier ``RRM``.

The first spin-labelling sitescan is performed with spin label ``mtsl`` on entity ``T4L`` and stored in file ``T4L_sites.lst``. The extension ``.lst`` is appended by default.
The demo then performs two spin-labeling sitescans on the RRM/stemloop complex. First, a scan with label ``mtsl`` is  performed on entity ``RRM``, generating site list ``RRM_sites``. The RRM is scanned only for substituting isoleucine residues by a spin label.
Second, a scan with label iodoacetamido-proxyl as a thiouracil label is performed on entity ``RRM``, generating site list ``SL_sites``. Only uracil nucleotides are considered for labeling.

The demo then generates two site pair lists, one from ``T4L_sites.lst`` referring to entityt ``T4L`` and stores it in file ``T4L_pairs.lst``. 
Site pairs are included only if the mean distance ranges between 15 and 80 Å. The second list combines labels from ``RRM_sites.lst`` with labels from ``SL_sites.lst``, referring to entity ``RRM``.
The pair list is stored as ``RRM_pairs.lst`` and again site pairs are included only if the mean distance ranges between 15 and 80 Å.

Then, basis file name ``T4L_distributions`` is declared for distance distribution output plots and file format ``pdf`` is selected for these plots.
Distributions are generated with label ``mtsl`` for entity ``T4L`` and are stored in numerical format with basis file name ``T4L_distribution``.
A list of 8 spin pairs is declared. Note that this computation of distance distributions does not require a previous site scan or pair list generation.

After this, a new basis file name ``RRM_distributions`` is declared for plots, which are again stored in ``pdf`` format. In this case, distributions are genarted for the existing pair list ``RRM_pairs``.
An orthogonal spin labelling scheme with labels ``mtsl`` and ``iap-4tu`` is used.

The final example of this demo imports PDB structure ``1BUY`` of erythropoietin into entity ``EPO``. At sitescan with ``mtsl`` is performed to provide list ``EPO_sites.lst``.
Site pairs are generated and scored from this list for elastic network modeling of conformation change of EPO. The list is stored as ``EPO_enm_pairs.lst``, considering only pairs with mean distance between 15 and 80 Å. 

Outputs
---------------------------------

Site and pair list are stored in a self-explaining human-readable format. Figures are stored as vector graphics in PDF format. 
Distance distributions are stored as comma-separated value (``.csv``) files with the distance in units of Å as the first column and probability density in units of :math:`\AA^{-1}` as the second column.