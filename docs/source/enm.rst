.. _enm:

ENM
==========================

This module generates ensembles from a known structure template of a protein or complex by deforming an elastic network model of the backbone to fit restraints.

``getpdb``
---------------------------------

Input of a raw ensemble (uniform populations) by reading a single PDB file. 

.. code-block:: matlab

    getpdb file

Arguments
    *   ``file`` - file name
Remarks
    *   the PDB file can contain several models (conformers) or a single one, by default, the first conformer is used
	
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
    *   '_m%i.pdb' is appended to the basis file name, where '%i' is the output model number
	
``conformer``
---------------------------------

Specifies input conformer from an input ensemble. 

.. code-block:: matlab

    conformer number

Arguments
    *   ``number`` - conformer number
Remarks
    *   if the key is missing, the first conformer is used
	
``remove``
---------------------------------

Removes a residue from the template. 

.. code-block:: matlab

    remove address

Arguments
    *   ``address`` - MMMx address of the residue to be removed
Remarks
    *   this can be used for removing cofactors of a template for fitting an apo structure
    *   use several ``remove`` keys if you wish to remove more than one residue
    *   for more complex manipulation of the template, use the Prepare module	
	
``chains``
---------------------------------

Use only some of the chains from the input structure. 

.. code-block:: matlab

    chains chain_id_1 [chain_id_2 ...]

Arguments
    *   ``chain_id_1`` - chain identifier of a chain to the included. e.g. 'A'
    *   ``chain_id_2 ...`` - optional chain identifiers of additional chains to the included. e.g. 'B' 'C'
Remarks
    *   this can be used if the template is a complex and the target is one of its components 
    *   there should be only one ``chains`` key or none at all
    *   by default, all chains are included in the template
	
``drag``
---------------------------------

Definition of residues dragged along with the elastic network. This is a block key with `n` lines for `n` residues. 

.. code-block:: matlab

    drag
       'address'
       ...
    .drag

Arguments
    *   ``address`` - address of a residue to be dragged along, e.g. `(A)501`
Remarks
    *   by default, only peptide chains are converted to a `C\alpha` elastic network model and deformed
    *   dragged residues are subjected to the same rotation and translation as the closest `C\alpha` atom  
    *   use this for ions and other cofactors
    *   it is advisable to refine the models afterwards

``ensemble``
---------------------------------

Specifies size of the output ensemble

.. code-block:: matlab

    ensemble size [uncertainty]

Arguments
    *   ``size`` - number of models in the output ensemble, defaults to 100
    *   ``uncertainty`` - optional uncertainty threshold, multiplier to standard deviation, defaults to 3
Remarks
    *   '_m%i.pdb' is appended to the basis file name, where '%i' is the output model number
    *   default uncertainty assumes subsequent ensemble fitting and contraction, use a lower value, if this is not intended
	
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

	
