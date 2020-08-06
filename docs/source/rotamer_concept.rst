.. _rotamer_concept:

Rotamer libraries
====================

Concept
---------------------------------

As previously established in MMM, site-directed labeling is modelled in MMMx based on pre-computed rotamer libraries of the label side groups.
Such libraries can be generated with classical atomic force fields via molecular dynamics simulations or Monte Carlo sampling of conformation space.
Briefly, the label side group is represented by a moderate number, typically between 100 and 10000, of rotameric states and their associated 
reference populations for a state that is non-interacting with the protein. Upon attachment to a labelling site, interaction with the protein is computed 
for a static (ensemble) structure of the latter assuming only non-bonded interaction via a Lennard-Jones potential. Flexibility of the protein as well as 
of the label beyond rotameric states is roughly modelled by an *f*-factor (`0 \le f \le 1`) that scales van-der-Waals radii of atoms in the underlying force field. 

Attachment results in rotation and translation of the rotamer coordinates and in reweighting of populations, assuming a Boltzmann distribution at 298 K. 
Distributions of properties can then be computed by population-weighted averaging. 

Implementation
--------------

The labeling concept is implemented by providing a set of rotamer libraries, a package for generating additional libraries for new labels, 
and by a function ``get_label`` for retrieving attributes of labels, such as distributions of label position and orientation.

Labelling is dynamic and virtual. Dynamic labeling means that a label *L* at residue *R* is computed at the time when the user of an entity first tries 
to retrieve attributes of this label at this particular site. The label attributes are then stored at residue level, whereas no atomic coordinates are 
generated at entity level (virtual labeling). Explicit coordinates of all atoms are generated only for clash tests, for saving label rotamers to a PDB file,
for visualizing the label in all-atom graphics, or when requested. Once computed, the explicit coordinates are stored in the atom coordinate array  of the
entity, while their indices are stored in the ``labels`` field of the labelled residue.

The following chemical modifications are implicit (d on need to be performed before labeling):

* mutation of any amino acid to cysteine or an unnatural amino acid for labeling

* conversion of uracil to 4-thiouracil

* conversion of phosphates to thiophosphates

.. _get_label:

The ``get_label`` function
----------------------------- 

Use this function if you want to label a residue with a particular label or want to retrieve attributes of a previously computed label.

.. code-block:: matlab

    argout = get_label(entity,label,attribute)
    [argout,exceptions] = get_label(entity,label,attribute)
    argout = get_label(entity,label,attribute,address)
    [argout,exceptions] = get_label(entity,label,attribute,address)


Parameters
    *   ``entity`` - entity in MMMx:atomic format
    *   ``label`` - label name
    *   ``attribute`` - see table below (string)
    *   ``address`` - MMMx residue address
Returns
    *   ``argout`` - output arguments (*M*-element cell array)
    *   ``exceptions`` - error message (1-element cell array)
	
Either a rotamer library for ``label`` must be registered with MMMx (see below) or the label name must be ``atom.<atname>``, where ``<atname>`` is an existing atom name for this residue.
	
**Attributes**
	
====================== =============================================== =================
Variable               Explanation                                     Type   
====================== =============================================== =================
``info``               rotamer library information                     struct
``info.tlc``           MMMx internal three-letter code                 string
``info.class``         label class, for instance ``nitroxide``         string
``info.smiles``        SMILES string representing label structure      string
``info.f_factor``      f-factor definition for attachment              double
``info.ref_pop``       reference populations for *C* rotamers          (*C*,1) double
``dihedrals``          *D* sidechain dihedrals for *C* rotamers        (*C,D*) double
``orientations``       Euler angles of molecular frame                 (*C*,3) double
``populations``        populations for *C* rotamers                    (*C*,1) double
``positions``          positions for *C* rotamers                      (*C*,3) double
====================== =============================================== =================
 
Note that SMILES strings for nitroxides tend to be interpreted as the corresponding hydroxylamines by some programs, notably by ChemDraw.
For libraries that contain several stereoisomers, the SMILES string refers to only one of them.

Set of libraries
-----------------

For the following labels, rotamer libraries are provided with MMMx:

=======  ===============================  ==============  ================ =========
Label    Synonyms                         Class           Attachment       Rotamers
=======  ===============================  ==============  ================ =========
``R1A``  ``mtsl``, ``mtssl``              nitroxide       cysteine         216
``R7A``  ``br-mtsl``, ``br-mtssl``        nitroxide       cysteine         216
``V1A``  ``v1``                           nitroxide       cysteine         72
``IA1``  ``ia-proxyl``                    nitroxide       cysteine         108
``MA1``  ``ma-proxyl``                    nitroxide       cysteine         108
``DZD``  ``dzd``                          nitroxide       cysteine         216
``DZC``  ``dzc``                          nitroxide       cysteine         216
``GDI``  ``iag``                          nitroxide       cysteine         2461
``GDM``  ``mag``                          nitroxide       cysteine         1367
``TUP``  ``iap-4tu``, ``iap-thiouracil``  nitroxide       4-thiouracil     72                          
``TUM``  ``mts-4tu``, ``mts-thiouracil``  nitroxide       4-thiouracil     192
``RTT``  ``r5-tpt``                       nitroxide       5'-thiophosphate 576
``R5P``  ``r5p``                          nitroxide       5'-thiophosphate 2048
``R3P``  ``r3p``                          nitroxide       3'-thiophosphate 512
``K1H``  ``HF-K1``                        nitroxide       unnatural aa     288
``NC1``  ``cNox@Tyr``                     nitroxide       tyrosine         128
``NX1``  ``lNox@Tyr``                     nitroxide       tyrosine         256   
``CNR``  ``CNC-NO``                       nitroxide       cofactor         144
``GMO``  ``dota-gd``                      gadolinium      cysteine         648
``GTO``  ``dtpa-gd``                      gadolinium      cysteine         2430
``M8D``  ``m8-dota-gd``                   gadolinium      cysteine         1944
``GPM``  ``gpymi-MTA``                    gadolinium      cysteine         432
``TMT``  ``tormyshev-trityl``             trityl          cysteine         3888
``HCU``  ``dHis-Cu``                      histidine       any amino acid   12
=======  ===============================  ==============  ================ =========

Label names (three-letter codes) and synonyms are case-insensitive. 
Note that gadolinium labels are sufficiently good approximations for other lanthanide labels with the same ligand,
for instance, for pseudo-contact shift (PCS) and paramagnetic relaxation enhncement (PRE) computations. 

-----------------------------

Atoms
---------

Use this function if you want to operate on atoms of all rotamers or on all atom locations.
Selections above atom level are ignored.

Use Matlab built-in function ``cell2mat`` for reforming output for B factor, charge, atomic number, and population into vectors. 
Do this only if you do not want to reassign them later after modification. If you do want to reassign, operate on the cell vectors. 
Note that MMMx supports only one B factor per atom, not distinct B factors for locations.

.. code-block:: matlab

    argout = get_atom(entity,attribute)
    [argout,exceptions] = get_atom(entity,attribute)
    argout = get_atom(entity,attribute,address)
    [argout,exceptions] = get_atom(entity,attribute,address)


Parameters
    *   ``entity`` - entity in MMMx:atomic format
    *   ``attribute`` - see table below (string)
    *   ``address`` - MMMx address for object selection, 'selected' or empty uses current selection
Returns
    *   ``argout`` - output arguments (*M*-element cell array)
    *   ``exceptions`` - error message, if attribute is not supported  (1-element cell array)
	
**Attributes**
	
====================== =============================================== ================================
Variable               Explanation                                     Type   
====================== =============================================== ================================
``bfactor``            crystallographic B factor, zero if unspecified  double
``charge``             atom charge, usually unspecified (zero)         int
``element``            atomic number                                   int8        
``info``               object information                              struct
``info.name``          atom name                                       string
``info.indices``       index vector (MMMx:atomic)                      (1,5) uint16 array
``info.atom_index``    index into atom array                           int
``population``         rotamer population or atom occupancy            double
``xyz``                Cartesian coordinates per location              (1,3) double
====================== =============================================== ================================ 

-----------------------------
   
Locations
---------

Use this function if you want to operate on selected rotamers or atom locations.
If the selection is on atom level and no rotamers are selected, only the first location or rotamer is referred to.
Selections above atom level are ignored.

.. code-block:: matlab

    argout = get_location(entity,attribute)
    [argout,exceptions] = get_location(entity,attribute)
    argout = get_location(entity,attribute,address)
    [argout,exceptions] = get_location(entity,attribute,address)


Parameters
    *   ``entity`` - entity in MMMx:atomic format
    *   ``attribute`` - see table below (string)
    *   ``address`` - MMMx address for object selection, 'selected' or empty uses current selection
Returns
    *   ``argout`` - output arguments (*M*-element cell array for *M* selected locations)
    *   ``exceptions`` - error message, if attribute is not supported  (1-element cell array)
	
**Attributes**
	
====================== =============================================== ================================
Variable               Explanation                                     Type   
====================== =============================================== ================================
``element``            atomic number                                   int8        
``info``               object information                              struct
``info.tag``           location tag, R# for a rotamer                  string, # is rotamer number
``info.indices``       index vector (MMMx:atomic)                      (1,5) uint16 array
``info.atom_index``    index into atom array                           int
``population``         rotamer population or atom occupancy            double
``xyz``                Cartesian coordinates per location              (1,3) double
====================== =============================================== ================================ 



