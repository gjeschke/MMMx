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

Computation of distance distributions
-------------------------------------

MMMx can compute distance distributions between two atoms, between an atom and a label, 
or between two labels for individual structural models or ensemble models. Thw two labels can be the same or they can differ.

Distance distributions for single conformers or ensembles given as PDB files can be computed, displayed, and saved 
by module :ref:`ExperimentDesign <experiment_design>`. They are automatically computed by module :ref:`EnsembleFit <ensemble_fit>`
for those site pairs where distance distribution restraints were specified. If you want to compute distance distributions between site pairs
that were not used in fitting, then first perform the ensemble fit and then run module :ref:`EnsembleFit <ensemble_fit>` again
with the fitted ensemble as input, the additional site pairs specified as distance distribution restraints, and the ``nofit`` keyword.