.. _experiment_design:

ExperimentDesign
==========================

This module helps in selecting site pairs for spectroscopic labeling. The sequence of command lines is relevant.
Features are demonstrated in examples ``demo_ExperimentDesign.mcx`` and ``demo_RBreference.mcx``. 
The latter example is a two-module pipeline that uses module :ref:`Prepare <prepare>` before ``EperimentDesign``. 

The following keywords are supported:


``import``
---------------------------------

Input of a raw ensemble (uniform populations) by reading a single PDB file. 

.. code-block:: matlab

    import file identifier

Arguments
    *   ``file`` - file name
    *   ``identifier`` - module-internal entity identifier of this ensemble
Remarks
    *   the PDB file can contain several models (conformers) or a single one
    *   for MMMx ensemble PDB files with population information in ``REMARK 400``, such information is discarded
	
``addpdb``
---------------------------------

Input of a single conformer or of a raw ensemble (uniform populations) by using wildcards in the file name. 

.. code-block:: matlab

    addpdb files identifier

Arguments
    *   ``files`` - file name, which can include wildcards (``*``, ``?``) for ensemble input
    *   ``identifier`` - module-internal entity identifier of this ensemble
Remarks
    *   for ensemble input with wildcards, all PDB files must contain conformers with the same primary structure
    *   a single input PDB file can contain several models (conformers), for clarity, rather use ``import`` for this case
    *   for MMMx ensemble PDB files with population information in ``REMARK 400``, such information is discarded

``getens``
---------------------------------

Input of an ensemble (with populations) from an MMMx :ref:`ensemble list <ensemble_list>` 

.. code-block:: matlab

    getens file identifier

Arguments
    *   ``file`` - file name, extension should be `.ens`
    *   ``identifier`` - module-internal entity identifier of this ensemble
Remarks
    *   keyword ``input`` is synonymous with ``getens``

``expand``
---------------------------------

Expand a rigid-body arrangement ensemble computed by the :ref:`Rigi <rigi>` module

.. code-block:: matlab

    expand file identifier

Arguments
    *   ``file`` - file name, extension should be `.mat`
    *   ``identifier`` - module-internal entity identifier of this ensemble
Remarks
    *   the whole ensemble will be built in memory, be cautious with very large ensembles

``plot``
---------------------------------

Requests that any generated plots are saved as graphics files.

.. code-block:: matlab

    plot file extension

Arguments
    *   ``file`` - basis file name, from which all plot file names are derived
    *   ``format`` - graphics format
Remarks
    *   possible graphics formats are 'pdb', 'png', 'epsc' (encapsulated postscript), 'jpg', 'bmp', 'emf' (enhanced metafile), 'tif'
    *   from experience, 'pdf' or 'epsc' is recommended for vector graphics and 'png' or 'tif' for bitmaps	
    *   if this command is missing, plots are not saved

``sitescan``, ``sitescan!``
---------------------------------

Spectroscopic-labeling site scans.

.. code-block:: matlab

    sitescan label entity outname [restypes [minrotamers [minpartf [chains]]]]

Arguments
    *   ``label`` - label, see :ref:`Rotamer libraries <rotamer_concept>` for list of available labels 
    *   ``entity`` - entity identifier specified in any of the input commands
    *   ``outname`` - file name for the site scan list output, extension is ``.lst``
    *   ``restypes`` - string of single-letter identifiers of residues to be considered, optional, defaults to 'CILMSTV'
    *   ``minrotamers`` - minimum number of rotamers for a site to be considered, optional defaults to 1
    *   ``minpartf`` - minimum partition function, optional, defaults to 0.1
    *   ``chains`` - restrict site scan to certain chains, string, such as 'AC', optional, defaults to '*' (all chains)
Remarks
    *   argument order matters, use defaults for earlier arguments, if you wish to deviate from a default in a later argument
    *   command ``sitescan`` considers only the first conformer in an ensemble, use ``sitescan!`` to scan all conformers

``pairlist``, ``pairlist!``
---------------------------------

Find feasible site pairs from site scan lists.

.. code-block:: matlab

    pairlist sitescan entity outname [rmin [rmax]]

Arguments
    *   ``sitescan`` - file name of a site scan list, generated, e.g., by the ``sitescan`` keyword, extension ``.lst`` 
    *   ``entity`` - entity identifier specified in any of the input commands
    *   ``outname`` - file name for the site scan list output, extension is ``.lst``
    *   ``rmin`` - minimum mean distance (Å) for a site pair to be considered, optional, defaults to 20 Å 
    *   ``rmax`` - maximum mean distance (Å) for a site pair to be considered, optional, defaults to 60 Å
Remarks
    *   argument order matters, use defaults for earlier arguments, if you wish to deviate from a default in a later argument
    *   command ``pairlist`` considers only pairs within the same conformer of an ensemble, use ``pairlist!`` to include inter-conformer pairs in computation of the distribution
    *   the label type is taken from the site scan list
    *   use ``hetpairlist`` or ``hetpairlist!`` if you want to combine sites with different labels

``hetpairlist``
---------------------------------

Find feasible site pairs from two site scan lists obtained for different labels (spectroscopically orthogonal labeling).

.. code-block:: matlab

    hetpairlist sitescan_1 sitescan_2 entity outname [rmin [rmax]]

Arguments
    *   ``sitescan_1`` - file name of the first site scan list, generated, e.g., by the ``sitescan`` keyword, extension ``.lst`` 
    *   ``sitescan_2`` - file name of the second site scan list 
    *   ``entity`` - entity identifier specified in any of the input commands
    *   ``outname`` - file name for the site scan list output, extension is ``.lst``
    *   ``rmin`` - minimum mean distance (Å) for a site pair to be considered, optional, defaults to 20 Å 
    *   ``rmax`` - maximum mean distance (Å) for a site pair to be considered, optional, defaults to 60 Å
Remarks
    *   argument order matters, use defaults for earlier arguments, if you wish to deviate from a default in a later argument
    *   command ``hetpairlist`` considers only pairs within the same conformer of an ensemble, use ``hetpairlist!`` to include inter-conformer pairs in computation of the distribution
    *   the label types are taken from the site scan lists

.. code-block:: matlab

    hetpairlist sitescan_1 sitescan_2 entity outname [rmin [rmax]]

``distributions``
---------------------------------

Compute, save, and plot distance distributions for directly specified site pairs or pair lists. This is a block key.

.. code-block:: matlab

    distributions labels entity outname [rmin [rmax [resolution]]]
       'address_1' 'address_2'
       []
       'pairlist'
    .distributions	

Arguments
    *   ``labels`` - label or labels for directly specified site pairs, see :ref:`Rotamer libraries <rotamer_concept>` for list of available labels 
    *   ``entity`` - entity identifier specified in any of the input commands
    *   ``outname`` - basis file name for the output distributions and plots
    *   ``rmin`` - minimum of distance axis, defaults to 10 Angstrom
    *   ``rmax`` - maximum of distance axis, defaults to 150 Angstrom
    *   ``resolution`` - resolution of distance axis, defaults to 0.5 Angstrom
    *   ``address_1`` - MMMx :ref:`residue address <MMMx_addresses>` for first site
    *   ``address_2`` - MMMx :ref:`residue address <MMMx_addresses>` for second site
    *   ``pairlist`` - file name of a site pair list obtained by commands ``pairlist`` or ``hetpairlist``
Remarks
    *   in the block you can specify as many site pairs and pair lists as you wish
    *   labels specified in the pair list(s) prevail over labels specified in the command, if these are inconsistent, it is reported in the log file
    *   if the same label is used at both sites, just provide this one; for orthogonal labeling, use the syntax ``label_1|label_2``
    *   if residue addresses do not contain a conformer specification, the distribution corresponds to all conformers

``trivariate``
---------------------------------

Compute and save trivariate distance distributions and plot their 2D and 1D sum projections for site triples. This is a block key.

.. code-block:: matlab

    trivariate labels entity outname
       'address_1' 'address_2' 'address_3'
       []
    .trivariate	

Arguments
    *   ``labels`` - label or labels for specified site triples, see :ref:`Rotamer libraries <rotamer_concept>` for list of available labels 
    *   ``entity`` - entity identifier specified in any of the input commands
    *   ``outname`` - basis file name for the output distributions and plots
    *   ``address_1`` - MMMx :ref:`residue address <MMMx_addresses>` for first site
    *   ``address_2`` - MMMx :ref:`residue address <MMMx_addresses>` for second site
    *   ``address_3`` - MMMx :ref:`residue address <MMMx_addresses>` for third site
Remarks
    *   in the block you can specify as many site triples as you wish
    *   if the same label is used at both sites, just provide this one; for orthogonal labeling, use the syntax ``label_1|label_2|label_3``
    *   output is large and thus in a binary Matlab file, variables are 'trivariate' for the 3D array and r_axis_1, 'r_axis_2', 'r_axis_3' for the three distance axes
    *   sequence of dimensions is 'site1-site2', 'site1-site3', 'site2-site3'
    *   use 'scipy.io.loadmat' from the SciPy library for importing to Python

``RigiFlex``
---------------------------------

Prepare RigiFlex rigid-body file and MMMx script template from AlphaFold prediction

.. code-block:: matlab

    RigiFlex UniProtID

Arguments
    *   ``UniProtID`` - sequence ID from UniProt, for which a prediction exists in the AlphaFold Protein Structure Database
Remarks
    *   output file names are automatically generated from the UniProt ID 
    *   if no prediction exists in the database, no output is generated and an error message is written to the log file
    *   the task also prepares spin-labelling site scan lists for all detected folded domains

``ENMpairs``
---------------------------------

Score site pairs for elastic network modeling.

.. code-block:: matlab

    ENMpairs sitescan entity outname [rmin [rmax]]

Arguments
    *   ``sitescan`` - file name of a site scan list, generated, e.g., by the ``sitescan`` keyword, extension ``.lst`` 
    *   ``entity`` - entity identifier specified in any of the input commands
    *   ``outname`` - file name for the site scan list output, extension is ``.lst``
    *   ``rmin`` - minimum mean distance (Å) for a site pair to be considered, optional, defaults to 20 Å 
    *   ``rmax`` - maximum mean distance (Å) for a site pair to be considered, optional, defaults to 60 Å
Remarks
    *   the output pair list is ordered by a score that predicts sensitivity of the pair to motion along the normal modes of the elastic network model 
    *   the input entity should be a single conformer, for an ensemble, only the first conformer is considered
    *   the label type is taken from the site scan list

``RBreference``
---------------------------------

Finds optimal reference sites in a rigid body by maximizing the area of the triangle spanned by three sites.

.. code-block:: matlab

    RBreference entity rmin rmax sitescan_1 [sitescan_2 ...]

Arguments
    *   ``entity`` - entity identifier specified in any of the input commands
    *   ``rmin`` - minimum mean distance (Å) for a site pair to be considered
    *   ``rmax`` - maximum mean distance (Å) for a site pair to be considered
    *   ``sitescan_1`` - file name of a site scan list, generated, e.g., by the ``sitescan`` keyword, extension ``.lst`` 
    *   ``sitescan_2 ...`` - optional file name(s) of further site scan lists, specify as many as you wish 
Remarks
    *   the output is provided in the log file in a format that can be used as input for the :ref:`Rigi <rigi>` module
    *   possibly, further chains contributing to the rigid body need to be added by the user
    *   see ``demo_RBreference.mcx`` for an example
