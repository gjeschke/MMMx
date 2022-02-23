.. _flex_RNA:

FlexRNA
==========================

This module generates single-stranded RNA linkers.

FlexRNA can have input arguments:

.. code-block:: matlab

    !flexRNA [coverage [models [maxtime]]]

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

``addpdb``
---------------------------------

Input of template conformers from PDB files. 

.. code-block:: matlab

    addpdb file

Arguments
    *   ``file`` - file name, can contain wildcards
Remarks
    *   use wildcard '*' for part of the filename to process all conformers from a previous step in the pipeline 
    *   without any input, FlexRNA generates free RNA chains
    *   use this command for attaching single-stranded RNA chains or linkers to existing structures
    *   in pipelines, use this command after previous Flex or FlexRNA modules
	
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
    *   use this command for processing of Rigi results by FlexRNA 
	
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
	
``sequence``
---------------------------------

nucleotide sequence for the RNA chain 

.. code-block:: matlab

    sequence nt_start nt_end seq

Arguments
    *   ``nt_start`` - number of the starting nucleotide, such as '4' 
    *   ``nt_end`` - number of the end residue, such as '8' 
    *   ``seq`` - sequence in single-letter format, such as 'UUCGA'
Remarks
    *   the sequence must consist of native nucleotides

``anchor_5p``
---------------------------------

5'-terminal anchor residue for the peptide chain 

.. code-block:: matlab

    anchor_5p  address

Arguments
    *   ``address`` - MMMx residue address, such as '(B)3' 
Remarks
    *   the addressed nucleotide must exist in the input conformers and must be a native nucleotide
    *   in pipelines with consecutive FlexRNA modules, address is affected by automatic chain identifier changes when chains are concatenated by linkers

``anchor_3p``
---------------------------------

3'-terminal anchor residue for the peptide chain 

.. code-block:: matlab

    anchor_3p address

Arguments
    *   ``address`` - MMMx residue address, such as '(C)9' 
Remarks
    *   the addressed nucleotide must exist in the input conformers and must be a native nucleotide
    *   in pipelines with consecutive FlexT=RNA modules, address is affected by automatic chain identifier changes when chains are concatenated by linkers

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
    *   if a residue is in the newly generated RNA, use only the residue number as its address
    *   the file name is optional, full distributions can be used
    *   if a full distribution is provided, ``rmean`` and ``rstd`` can be skipped
    *   distance distribution restraints are always treated as full distribution, if only ``rmean`` and ``rstd`` are provided, the distance is computed
	*   test of distance distribution restraints is done with full models and based on the overlap metric
 
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
