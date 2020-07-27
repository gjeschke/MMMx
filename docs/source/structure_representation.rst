Structure representation
=========================

Protein database (PDB) format
-----------------------------

Structural biology has inherited structure representation from protein crystallography with one addition from NMR structure determination.

The main aspects of this representation are:

1) Three-dimensional structure is defined by Cartesian atom coordinates in units of Angstroem

2) Uncertainty about structure is quantified by isotropic or anisotropic distribution of atom positions (B factors) for crystal structures.
   For NMR structures, it is represented by presenting several models (typically 20).
   
3) Genuine disorder is represented by alternate atom locations with associated populations or missing atoms in crystal structures.
   In NMR structures, genuine disorder contributes to structure variation between models and is generally indistinguishable from uncertainty.

In the end, this boils down to :math:`n_{\mathrm{coor},k}` coordinate specifications :math:`(x_k,y_k,z_k)` for atom :math:`k` of the entity, where

.. math::
   n_{\mathrm{coor},k} = \sum_{m = 1}^M \sum_{l = 1}^{L_{k,m}} 1 

where :math:`M` is the number of models and :math:`L_{k,m}` is the number of alternate locations of atom :math:`k` in model :math:`m`.

Whereas the alternate locations have associated populations :math:`p_{l,m}` with :math:`\sum_{l = 1}^{L_{k,m}} p_{l,m} = 1`, all models are assumed to be equally populated.

This format is supported by virtually all software packages for structural biology and will not vanish. MMMx can read and write it. 
It is not a suitable representation for the structure of entities with functionally relevant disorder.

PDBx/mmCIF format
------------------------------------

This format `has formally replaced the PDB format at protein database in 2014`__. It allows for much more annotation, is nicely extensible, and less strict on line formatting.

.. __: http://mmcif.wwpdb.org/docs/faqs/pdbx-mmcif-faq-general.html 

It does not differ from the PDB format with respect to structure representation.

MMMx internal representations
-----------------------------

As a multi-scale modeller, MMMx needs extensible and varied representations of entity structure. Two representations exist in the first release:

1) :ref:`MMMx|atomic<MMMx_atomic>`: This representation is closest to PDB representation, but allows for non-uniform conformer populations and sidegroup rotamers with associated populations. 

2) :ref:`MMMx|RigiFlex<MMMx_RigiFlex>`: This representation keeps full information on rigid-body arrangement for models where some chains or chain segments are specified as rigid bodies. 
   It is also more memory economic for such models, which can matter for large ensembles.
   
All MMMx internal representations can be converted to the PDB representation, although with some loss of information or large inflation of file size and memory requirements.

Conversion of representations can be sequential. 
For instance, the ``MMMx|RigiFlex`` representation is first expanded to the ``MMMx|atomic`` representation and then projected onto the PDB representation:

.. math::
   \mathrm{MMMx}|\mathrm{RigiFlex} \rightarrow \mathrm{MMMx}|\mathrm{atomic} \rightarrow \mathrm{PDB}
   
In this sequence, ``RigiFlex`` presents the highest level of description (most information) and ``PDB`` the lowest. Backwards conversion is not usually possible.
Internally, MMMx represents structure on the highest possible level.

It follows that ChimeraX holds the master copy of an entity structure only if the structure originally came from a PDB file or can be represented in PDB without loss of information.
Otherwise, MMMx holds the master copy.



