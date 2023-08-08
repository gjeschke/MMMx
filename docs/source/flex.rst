.. _flex:

Flex
==========================

This module generates peptide linkers with partial or no secondary structure.

Flex can have input arguments:

.. code-block:: matlab

    !flex [coverage [models [maxtime]]]

Arguments
    *   ``coverage`` - controls fraction of conformation space covered by the ensemble, defaults to 0.5  
    *   ``models`` - number of models per input conformer, defaults to 1  
    *   ``maxtime`` - maximum time spent on one model, defaults to 1 h 
Remarks
    *   coverage is not rigorously defined, but only an empirical control parameter
    *   the larger coverage is, the broader is the ensemble	
    *   change default of coverage only if you have very good reasons
    *   note that total maximum time can be as long as the product of `maxtime` abd the number of input conformers 

The following keywords are supported:

``a_prop``
---------------------------------

Definition of `\alpha`-helix propensities. This is a block key with `n` lines for `n` restraints. 

.. code-block:: matlab

    a_prop
       'resnum' 'propensity'
       ...
    .a_prop

Arguments
    *   ``resnum`` - residue number, e.g. `120`
    *   ``propensity`` - `\alpha`-helical propensity, number between 0 (none) and 1 (always `\alpha`-helical)
Remarks
    *   secondary structure propensities can be specified only for residues in the newly generated peptide
    *   `\alpha`-helical propensity is realized by imposing backbone dihedrals `\psi` and `\phi` corresponding to `\alpha`-helical structure

``acceptance``
---------------------------------

Controls acceptance threshold. 

.. code-block:: matlab

    acceptance threshold [mode]

Arguments
    *   ``threshold`` - fraction of models that should be accepted
    *   ``mode`` - acceptance mode, can be 'uniform' or 'individual', default is 'uniform'  
Remarks
    *   this option requires that full distance distributions are used as restraints
    *   with the acceptance key switches, a variant of von-Neumann rejection sampling is used
    *   by default, Gaussian restraints are used even if full distributions are provided
    *   the higher `threshold` is, the faster is model generaton, but the worse is agreement of the raw ensemble with distance distributions
    *   use higher `threshold` if model yield is too low for the downstream part of the pipeline 

``addpdb``
---------------------------------

Input of template conformers from PDB files. 

.. code-block:: matlab

    addpdb file

Arguments
    *   ``file`` - file_name, can contain wildcards
Remarks
    *   use wildcard '*' for part of the filename to process all conformers from a previous step in the pipeline 
    *   without any input, Flex generates free peptide chains
    *   use this command for attaching flexible peptide chains or linkers to existing structures
    *   in pipelines, use this command after previous Flex or FlexRNA modules
	
``b_prop``
---------------------------------

Definition of `\beta`-strand propensities. This is a block key with `n` lines for `n` restraints. 

.. code-block:: matlab

    b_prop
       'resnum' 'propensity'
       ...
    .b_prop

Arguments
    *   ``resnum`` - residue number, e.g. `147`
    *   ``propensity`` - `\beta`-strand propensity, number between 0 (none) and 1 (always `\beta`-strand)
Remarks
    *   secondary structure propensities can be specified only for residues in the newly generated peptide
    *   `\beta`-strand propensity is realized by imposing backbone dihedrals `\psi` and `\phi` corresponding to `\beta`-strand structure

``c_anchor``
---------------------------------

C-terminal anchor residue for the peptide chain 

.. code-block:: matlab

    c_anchor address

Arguments
    *   ``address`` - MMMx residue address, such as '(D)121' 
Remarks
    *   the addressed residue must exist in the input conformers and must be a native amino acid
    *   in pipelines with consecutive Flex modules, address is affected by automatic chain identifier changes when chains are concatenated by linkers

``c_prop``
---------------------------------

Definition of cis-propensities. This is a block key with `n` lines for `n` restraints. 

.. code-block:: matlab

    c_prop
       'resnum' 'propensity'
       ...
    .c_prop

Arguments
    *   ``resnum`` - residue number, e.g. `211`
    *   ``propensity`` - cis-propensity, number between 0 (always trans) and 1 (always cis)
Remarks
    *   cis-propensities can be specified only for residues in the newly generated peptide
    *   cis-propensity is realized by imposing backbone dihedral `\omega = 0^\circ` corresponding to  a cis-residue
    *   cis conformation usually occurs only for proline residues

``clashtest``
---------------------------------

Number of generated residues after which intermediate clashtests are performed

.. code-block:: matlab

    clashtest spacing

Arguments
    *   ``spacing`` - spacing between intermediate clashtests during backbone generation, defaults to 10000 (practically never)
Remarks
    *   change default only if you suspect a problem that can be solved this way
    *   usually, intermediate clashtests slow down model generation

``ddr``
---------------------------------

Definition of distance distribution restraints. This is a block key with `n` lines for `n` restraints. 

.. code-block:: matlab

    ddr label_1 [label_2]
       'address_1' 'address_2' 'rmean' 'rstd' [@'fname']
       ...
    .ddr

Arguments
    *   ``label_1``, ``label_2`` - label types, e.g. `mtsl`, `dota-gd`
    *   ``address_1``, ``address_2`` addresses of the two labelled sites, e.g., `(A)16`, `107`
    *   ``rmean`` mean distance in Angstroem, e.g. `32.5`
    *   ``rstd`` standard deviation in Angstroem, e.g. `15.5`
    *   ``fname`` optional file name of the distance distribution 
Remarks
    *   if both labels are the same, it is sufficient to specify the label type once
    *   use separate 'ddr' blocks for each label combination
    *   if a residue is in the newly generated peptide, use only the residue number as its address
    *   the file name is optional, full distributions can be used
    *   if a full distribution is provided, ``rmean`` and ``rstd`` can be skipped, these parameters are then automatically computed from the distribution
    *   for monomodal distributions, the advantage of using full distributions in terms of ensemble quality is (at best) minor 
    *   using full distributions provides more convenient control over model yield with the 'acceptance' keyword
 
``depth``
---------------------------------

Definition of bilayer immersion depth restraints. This is a block key with `n` lines for `n` restraints. 

.. code-block:: matlab

    depth label
       'resnum' 'rmean' 'rstd'
       ...
    .depth

Arguments
    *   ``label`` - label types, e.g. `CA` for Calpha
    *   ``resnum`` - residue number of the site, e.g., `3`
    *   ``rmean`` mean distance from bilayer central plane in Angstroem, e.g. `20`
    *   ``rstd`` standard deviation if the distribution in Angstroem, e.g. `15.5`
    *   ``fname`` file name of the distance distribution 
Remarks
    *   input structures must be in a frame where the bilayer normal is the z axis, use Prepare
    *   use `CA` as label identifier if you are unsure
    *   use separate 'depth' blocks for different labels
    *   depth restraints can be specified only for sites in the newly generated peptide
    *   use a negative argument instead of `rmean` for specifying a lower bound
    *   use a negative argument instead of `rstd` for specifying an upper bound
	
``expand``
---------------------------------

Input and expansion of rigid-body arrangements. 

.. code-block:: matlab

    expand [file]

Arguments
    *   ``file`` - optional fle name for rigid-body arrangements
Remarks
    *   without input argument, the output of a previous Rigi module in the pipeline is expanded 
    *   input file format is the Matlab output format of Rigi
    *   use this command for processing of Rigi results by Flex 
	
``getpdb``
---------------------------------

Input of a raw ensemble (uniform populations) by reading a single PDB file. 

.. code-block:: matlab

    getpdb file

Arguments
    *   ``file`` - file name
Remarks
    *   the PDB file can contain several models (conformers) or a single one
    *   for MMMx ensemble PDB files with population information in ``REMARK 400``, such information is read
	
``loose``
---------------------------------

Switches off sidechain clash test

.. code-block:: matlab

    loose

Remarks
    *   this option is intended only for cases where model generation is extremely slow or impossible otherwise
    *   do not use models obtained with the `loose` option without subsequent refinement (e.g. using YasaraRefine)
    *   models may clash so strongly that refinement with other programs fails

``n_anchor``
---------------------------------

N-terminal anchor residue for the peptide chain 

.. code-block:: matlab

    n_anchor address

Arguments
    *   ``address`` - MMMx residue address, such as '(A)89' 
Remarks
    *   the addressed residue must exist in the input conformers and must be a native amino acid
    *   in pipelines with consecutive Flex modules, address is affected by automatic chain identifier changes when chains are concatenated by linkers

``oligomer``
---------------------------------

Definition of oligomer distance distribution restraints. This is a block key with `n` lines for `n` restraints. 

.. code-block:: matlab

    oligomer label n
       'resnum' 'rmean' 'rstd' [@'fname']
       ...
    .oligomer

Arguments
    *   ``label`` - label types, e.g. `ia-proxyl`
    *   ``n`` - number of symmetry-related protomers in the oligomer, e.g. `3`
    *   ``resnum`` - residue number of the site, e.g. `7`
    *   ``rmean`` mean distance in Angstroem, e.g. `32.5`
    *   ``rstd`` standard deviation in Angstroem, e.g. `15.5`
    *   ``fname`` file name of the distance distribution 
Remarks
    *   input structures must be in a frame where the Cn symmetry axis is the z axis, use Prepare
    *   use separate 'oligomer' blocks for different labels
    *   oligomer restraints can be specified only for sites in the newly generated peptide
    *   the file name is optional, full distributions can be used
    *   the use of full distributions is implemented, but has not yet been tested in detail
	
``parallel``
---------------------------------

Controls parallelization of conformer generation. 

.. code-block:: matlab

    parallel trials

Arguments
    *   ``trials`` - number of trials computed in parallel before analysis, defaults to 100
Remarks
    *   change default only if you have a very good reason

``save``
---------------------------------

Specifies basis name for saving output conformers 

.. code-block:: matlab

    save file [[pdb_id] chain_id]

Arguments
    *   ``file`` - basis file name 
    *   ``pdb_id`` - optional four-letter (pseudo) PDB identifier
    *   ``chain_id`` - optional chain identifier
Remarks
    *   '_i%i_m%i.pdb' is appended to the basis file name, the first '%i' is input conformer number, the second '%i' is the model number for this input
    *   if a chain identifier is provided, a free-standing peptide gets this identifier
	
``scwrl4``
---------------------------------

Specify full SCWRL4 pathname. 

.. code-block:: matlab

    scwrl4 pathname

Arguments
    *   ``pathname`` - full path to the SCWRL4 executable
Remarks
    *   this is required on Linux systems, where Matlab does not find SCWRL4 even if it is on the Matlab path 
	
``sequence``
---------------------------------

amino acid sequence for the peptide chain 

.. code-block:: matlab

    sequence res_start res_end seq

Arguments
    *   ``res_start`` - number of the starting residue, such as '90' 
    *   ``res_end`` - number of the end residue, such as '120' 
    *   ``seq`` - sequence in single-letter format, such as 'RSGRGTGRGGGGGGGGGAPRGRYGPPSRRSE'
Remarks
    *   the sequence must consist of native amino acids

``skipto``
---------------------------------

Skips input conformers. 

.. code-block:: matlab

    skipto first

Arguments
    *   ``first`` - first input conformer for which models are generated
Remarks
    *   by default, there is no skipping
    *   this can be used after a crash or job timeout

``verbose``
---------------------------------

Sets verbose mode. 

.. code-block:: matlab

    verbose [trials]

Arguments
    *   ``cycles`` - number of Monte carlo trials after which new verbose information is written to logfile
Remarks
    *   by default, verbose is off
    *   verbose without argument has a default of 200 trials
    *   verbose writes time per generated model, an estimate of remaining computation time, and statistics on the reasons for failed trials
