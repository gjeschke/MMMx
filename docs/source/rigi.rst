.. _rigi:

Rigi
==========================

This module generates ensembles of rigid-body arrangements. Sequence of commands is irrelevant.

Features are demonstrated in example ``demo_Rigi.mcx``, which uses module ``Prepare`` to generate a rigid-body
template from a structure downloaded from the PDB server. Bare templates can also be prepared from AlphaFold predictions downloaded
from the AlphaFold Protein Structure Database by UniProt identifier with the ``Rigi`` keyword in ``ExperimentDesign`` 

The following keywords are supported:

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
    *   ``address_1``, ``address_2`` addresses of the two labelled sites, e.g., `(A)16`, `(D)148`
    *   ``rmean`` mean distance in Angstroem, e.g. `32.5`
    *   ``rstd`` standard deviation in Angstroem, e.g. `15.5`
    *   ``fname`` file name of the distance distribution 
Remarks
    *   if both labels are the same, it is sufficient to specify the label type once
    *   use separate 'ddr' blocks for each label combination
    *   the file name is optional, full distributions are not currently used
	
``maxsize``
---------------------------------

Specifies the maximum size (extension) of a rigid-body arrangement  

.. code-block:: matlab

    maxsize size

Arguments
    *   ``size`` - maximum distance between reference points in a rigid-body arrangement, defaults to 180 Angstroem
  
``maxtime``
---------------------------------

Maximum run time in hours.  

.. code-block:: matlab

    maxtime t

Arguments
    *   ``t`` - approximate maximum run time allowed for Rigi
Remarks
    *   computation can take somewhat longer, because runtime ist tested only execution of full parallel blocks (10,000 trials)
    *   if Rigi is interrupted by timeout, sampling is not exhaustive
    *   for production runs, reduce maxtrials until you can complete the exhaustive sampling in a reasonable time
    *   default is 48 h

``maxtrials``
---------------------------------

Maximum number of trials in exhaustive search.  

.. code-block:: matlab

    maxtrials T

Arguments
    *   ``T`` - maximum number of trials, determines, how long an exhaustive search takes
Remarks
    *   this must be specified, start with 10,000 if you are not sure and see how long it takes
    *   actual number of trials is a product of integer digital resolutions for all core restraints and may be smaller
    *   the larger ``T``, the better the resolution of the exhaustive search
    *   computation time is linear in actual number of trials
	
``models``
---------------------------------

Maximum number of rigid-body arrangements in output.  

.. code-block:: matlab

    models M

Arguments
    *   ``M`` - maximum number of models that are returned
Remarks
    *   if exhaustive sampling yields less models, this setting has no effect
    *   if exhaustive sampling yields more models, the solutions are hierarchically cluster to ``M`` models
    *   use this, if Rigi returns too many models for further processing
    *   the default is 20,000

``nlink``
---------------------------------

Nucleotide link. This is a block key with `n` lines for `n` links. 

.. code-block:: matlab

    nlink 
       'address_1' 'address_2' 'nucleotides' 'length'
       ...
    .nlink

Arguments
    *   ``address_1``, ``address_2`` addresses of the two anchor nucleotides, e.g., `(B)3`, `(C)6`
    *   ``nucleotides`` linker segments, number of missing nucleotides + 1
    *   ``length`` maximum length in Angstroem, up to 6*nucleotides
Remarks
    *   the anchor nucleotides must exist in rigid bodies
    *   slightly shorter lengths, e.g. 16 instead of 18 for two missing nucleotides, improve success rate of FlexRNA
	
``parallel``
---------------------------------

Specifies number of trials in a parallel block  

.. code-block:: matlab

    parallel ptrials

Arguments
    *   ``ptrials`` - number of trials performed in parallel, defaults to 10,000
Remarks
    *   change from default only if you have a very good reason   	
	
``plink``
---------------------------------

Peptide link. This is a block key with `n` lines for `n` links. 

.. code-block:: matlab

    plink 
       'address_1' 'address_2' 'residues' 'length'
       ...
    .plink

Arguments
    *   ``address_1``, ``address_2`` addresses of the two anchor residues, e.g., `(A)89`, `(D)121`
    *   ``residues`` linker segments, number of missing residues + 1
    *   ``length`` maximum length in Angstroem, up to 3.8*residues
Remarks
    *   the anchor residues must exist in rigid bodies
    *   slightly shorter lengths improve success rate of Flex

``probability``
---------------------------------

Specifies the probability covered by the RBA ensemble  

.. code-block:: matlab

    probability p

Arguments
    *   ``p`` - number between 0 and 1, defaults to 0.5
Remarks
    *   change from default only if you have a very good reason   	

``rbtemplate``
---------------------------------

Input of a rigid-body template file.  

.. code-block:: matlab

    rbtemplate file

Arguments
    *   ``file`` - PDB file name, must include extension `.pdb`, otherwise download from the PDB server is attempted
Remarks
    *   the command ``addpdb`` is synonymous with ``rbtemplate``
    *   you can use the ``merge`` command in module :ref:`Prepare <prepare>` to generate a rigid-body template file

``resolution``
---------------------------------

Sets a resolution limit for the exhaustive search of RBA arrangement space  

.. code-block:: matlab

    resolution res

Arguments
    *   ``res`` - resolution in Angstroem, defaults to 3 Angstroem
Remarks
    *   actual resolution can be larger, but not smaller   		    
    *   do not change from default, if you are not sure  		
	
``rigid``
---------------------------------

Definition of rigid bodies. This is a block key with `n` lines for `n` rigid bodies. 

.. code-block:: matlab

    rigid chain_1 [chain_2 ...]
       'address_1' 'label_1'
       'address_2' 'label_2'
       'address_3' 'label_3'
    .rigid

Arguments
    *   ``chain_1``, ``chain_2``, ... - chains belonging to this rigid body, example `(A)` `(B)`
    *   ``address_1`` address of the first reference point, e.g., `(A)16`
    *   ``label_1`` label type for the first reference point, e.g. `mtsl`
Remarks
    *   use exactly three reference points 
    *   at least one chain and as many chains as needed can belong to one rigid body
    *   reference point addresses mus be in one of the chains that belong to the rigid body
    *   you can define as many rigid bodies as you need, but computational effort increases exponentially
    *   there must be at least two rigid bodies
	
``save``
---------------------------------

Specifies name for saving the output in MMMx:rigid_body format.  

.. code-block:: matlab

    save fname

Arguments
    *   ``fname`` - file name for output, extension ``.mat`` is appended if none
Remarks
    *   if not present, output is automatically save to ``MMMx_rigi.mat``   

``savepdb``
---------------------------------

Specifies basis name for saving individual rigid-body arrangements to PDB files 

.. code-block:: matlab

    savepdb bname

Arguments
    *   ``bname`` - basis file name for PDB files, ``_rba_%i.pdb`` is appended, where ``%i`` denotes the number of the RBA
Remarks
    *   if not present, no individual PDB files are saved   	
	
``separate off``
---------------------------------

Turn off automatic separation of rigid bodies in the template file. 

Arguments
    *   none, ``separate off`` is the only syntax that has an influence
Remarks
    *   do this only if the rigid bodies are already well separated in the template
    *   if you wish to superimpose the ensemble onto a template, this is better done in module :ref:`EnsembleAnalysis<ensemble_analysis>`  

``superimpose``
---------------------------------

Superimpose all rigid-body arrangements at one rigid body.  

.. code-block:: matlab

    superimpose rigid_body

Arguments
    *   ``rigid_body`` - number of the rigid body at which arrangments are superimposed
Remarks
    *   the number corresponds to the sequence of ``rigid`` blocks in the control file   

``xl_percentage``
---------------------------------

Specifies the percentage of crosslink restraints that must be fulfilled for an RBA to be accepted  

.. code-block:: matlab

    xl_percentage p

Arguments
    *   ``p`` - number between 0 and 100, defaults to 30%
Remarks
    *   there is few experience what is appropriate for flexible systems   	
    *   directive has no effect, if no crosslinks are specified   	

``xlink``
---------------------------------

Crosslink. This is a block key with `n` lines for `n` links. 

.. code-block:: matlab

    xlink 
       'address_1' 'address_2' 'distance'
       ...
    .xlink

Arguments
    *   ``address_1``, ``address_2`` addresses of the two crosslinked residues
    *   ``distance`` maximum distance between CA atoms for the crosslink deemed to be possible
Remarks
    *   the crosslinked residues must exist in rigid bodies
    *   only a certain percentage of crosslinks needs to be fulfilled for the RBA to be valid, default is 30%
    *   use 'xlink_percentage' to set this percentage
    *   this feature is not well tested
    *   it is hard to predict which percentage of crosslinks should be fulfilled in any given arrangement
