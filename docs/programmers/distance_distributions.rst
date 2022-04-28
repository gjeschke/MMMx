.. _distance_distributions:

Distance distributions
==========================

Reasons why distances are distributed
-------------------------------------

In atomic-resolution structures obtained by x-ray crystallography, distances between atoms are defined with a 
precision given by the B factors. Parts of the macromolecules that are too disordered even in a crystal are not specified. 
Such disorder may arise from sidechains adopting several rotameric states or from coexistence of many different backbone conformations.

Atomic-resolution structures obtained by NMR are represented by several models, typically 20 models. For small globular proteins,
which are the mainstay of solution NMR, usually a core of the protein is defined with a atom position variance not much different from the
properly scaled B factor in a crystal structure of the same protein. Terminal loops, internal loops, and sometimes whole domains may
exhibit much more variation between models. Usually, this variation can only be interpreted in a binary way: this section is disordered.
The extent of disorder cannot be specified, as the contribution from a lack of restraints cannot safely be separated from the contribution due
to coexistence of several conformers.

Distance distributions, as they are accessible by pulse dipolar spectroscopy (PDS) EPR techniques, have the potential for at least
partial separation of the two contributions. In that, two problems arise from the necessity to introduce labels. First, the 
label might influence the ensemble of backbone conformers and, second, the non-native spin label side chain may have a broader distributions
of rotameric states than a native sidechain. The first problem appears to be minor compared to the second one at least in ordered cores.
More caution and more research may be required for disordered sections. It is generally advisable to integrate PDS data with data from other
techniques and to assess consistency between the individual restraint sets.

The second problem is addressed in MMMx by mdelling the conformational distribution of the spin label sidechain :ref:`with a rotamer library<rotamer_concept>`.
The same concept is applied to chromophore labels for Förster resonance energy transfer (FRET). 
Uncertainties in this modelling limit the resolution of structures that rely on PDS or FRET restraints and thus the minimal extent
of backbone disorder that can safely be recognized as backbone disorder. For the most widely applied methanthiosulfonate spin label (MTSL), 
the resolution limit is about 2-3 Å.

The attribute ``info`` cannot generally be modified, as this could compromise integrity of addressing or indexing the entity. 
There are exceptions for fields of info on some hierarchy levels.

Computation of distance distributions
-------------------------------------

MMMx can compute distance distributions between two atoms, between an atom and a label, 
or between two labels for individual structural models or ensembles models. All these computations use the same function:

.. code-block:: matlab

    [r_axis,distribution] = distance_distribution(entity,site1,label1,site2,label2)
    [r_axis,distribution,entity] = distance_distribution(entity,site1,label1,site2,label2)
    [r_axis,distribution,entity,exceptions] = distance_distribution(entity,site1,label1,site2,label2)
    [r_axis,distribution] = distance_distribution(entity,site1,label1,site2,label2,options)
    [r_axis,distribution,entity] = distance_distribution(entity,site1,label1,site2,label2,options)
    [r_axis,distribution,entity,exceptions] = distance_distribution(entity,site1,label1,site2,label2,options)

Parameters
    *   ``entity`` - entity in MMMx:atomic format
    *   ``site1`` - residue address of the first site
    *   ``label1`` - label name (see below) or ``atom.`` *atomname* for an atom in that residue (first site)
    *   ``site2`` - residue address of the second site
    *   ``label2`` - label name (see below) or ``atom.`` *atomname* for an atom in that residue (second site)
	*   ``options`` - computation options, see tble below
Returns
    *   ``r_axis`` - distance axis (Å)
    *   ``distribution`` - distance distributions
	*   ``entity`` - input entity augmented by newly computed labels
    *   ``exceptions`` - error messages  (cell array)
	
**Options**
	
====================== =============================================== ================================
Variable               Explanation                                     Type   
====================== =============================================== ================================
``rmin``               minimum distance, default: 10 Å                 double
``rmax``               maximum distance, default: 150 Å                double
``resolution``         resolution, default: 0.5 Å                      double
``units``              ``probability`` or  ``density`` (default)       string                    
``coupled``            for ensemble computations, ``coupled = true``   boolean
                       implies that both sites are in the same model
                       (default), else distributions between sites in
                       different models are also added
``smoothing``          Standard deviation for Gaussian smoothing,      double
                       defaults to twice ``options.resolution``
====================== =============================================== ================================ 

For ``options.units = 'probability'``, the sum of the distribution is unity.
For ``options.units = 'density'``, elements of the distribution vector have the unit 1/Å.  

It is good practice to receive the updated entity, as this saves time in later computations involving the same label(s).

The possible choices for ``label1`` and ``label2`` are given by the :ref:`implemented rotamer libraries<label_set>`. The syntax:

.. code-block:: matlab

    [r_axis,distribution] = distance_distribution(entity,'{*}(A)78','atom.CA','{*}(A)135','atom.CA');

returns the CA-CA distance distribution for residues 78 and 135 in chain A over all conformers (models) in the entity.
If only atoms are addressed, it is not necessary to receive the entity, as it remains unchanged.

The syntax:

.. code-block:: matlab

    [r_axis,distribution,entity] = distance_distribution(entity,'{3}(A)29','mtsl','{3}(B)507','dota-gd');

returns the distance distribution within conformer 3 between MTSL attached to residue 29 in chain A and a Gd(DOTA) label attached to residue 507 in chain B.

If the site address misses a conformer specification, conformer 1 is assumed. If the function fails, the distribution ist empty.

