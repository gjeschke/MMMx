.. _object_access:

Get object attributes
==========================

Concept
---------------------------------

Object access functions return or set object attributes, where an object is a conformer, a chain, a residue, a rotamer, an atom, or an atom location. 
They operate on an ``entity``, which is an MMMx ensemble structure representation in :ref:`MMMx|atomic<MMMx_atomic>` or :ref:`MMMx|RigiFlex<MMMx_RigiFlex>` format.

Wherever possible, methods should operate on an entity only through object access functions. In order to speed up access to coordinates, 
an additional function ``get_coor`` exists for fast retrieval of atom coordinates, atomic numbers, and corresponding indices into MMMx:atomic tables. 
This function works for objects selected on any level of structural hierachy and combinations of objects at various levels.
The pendant ``set_coor`` uses the indices returned by ``get_coor`` for reassignment of processed coordinates.

``get`` functions retrieve attributes, whereas ``set`` functions set or change attributes. They are specific to an object hierarchy level in order to avoid unintended behaviour.
The objects are selected by an :ref:`MMMx address<MMMx_addresses>`. The address can also be ``selected`` for accessing all currently selected objects on this hierarchy level. 

Generic syntax
--------------

.. code-block:: matlab

    [argout,exceptions] = get_"object"(entity,attribute,address)
	 
where the ``attribute`` strings are specific to the "object" hierarchy level 
("chain", "residue", "atom", "location", "conformer", see below) and argout is a cell vector, 
whose length is (most of the time) the number of objects of this hierarchy level that were selected by ``address``.
If the address is ``selected`` or empty, the ``get`` functions operate on objects currently selected in the entity.

Error messages or warnings are reported as MException objects in cell array ``exceptions``. 

Chain
---------

Use this function if you want to operate exclusively on chain level.
Selections below chain level are ignored.

.. code-block:: matlab

    argout = get_chain(entity,attribute)
    [argout,exceptions] = get_chain(entity,attribute)
    argout = get_chain(entity,attribute,address)
    [argout,exceptions] = get_chain(entity,attribute,address)


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
``info``               object information                              struct
``info.name``          chain tag                                       int
``info.type``          chain/molecule type, 'peptide', 'peptide+',     string
                       'DNA', 'DNA+', 'RNA', 'RNA+', 'HET', 'water'
``populations``        populations for all *C* conformers              (*C*,1) double
====================== =============================================== ================================ 
 
The ``+`` types indicate that a chain contains peptide, DNA, or RNA sections, but additionally other "residues".
Structures with such "chains" are frequent in the PDB. In native MMMx structures, other molecules or ions should 
be combined with a chain if and only if they are coordinated to the chain and coordination geometry is considered to
be immutable within the model approximation.   

Residue
---------

Use this function if you want to operate exclusively on residue level.
Selections above or below residue level are ignored.

Use Matlab built-in function ``cell2mat(argout)`` for reforming output for ``pointindices`` into a vector. 
To reformat the point coordinates ``pointcoor`` or the sheet information ``sheet`` into arrays, use

.. code-block:: matlab
    
    arrayout = argout(~cellfun('isempty',argout));
    outputs = length(arrayout);
    arrayout = cell2mat(arrayout);
    arrayout = reshape(arrayout,*n*,outputs).';
	
where ``*n* = 3`` for point coordinates and ``*n* = 2`` for sheet information. 
Do this only if you do not want to reassign these values after modification.
If you do want that, operate on the cell vectors.

Populations of rotamers have to be processed residue-by-residue, since the number of rotamers can differ between residues. 

.. code-block:: matlab

    argout = get_residue(entity,attribute)
    [argout,exceptions] = get_residue(entity,attribute)
    argout = get_residue(entity,attribute,address)
    [argout,exceptions] = get_residue(entity,attribute,address)


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
``dssp``               DSSP secondary structure assignment             char
``info``               object information                              struct
``info.number``        residue number                                  int
``info.tlc``           three-letter code/PDB residue tag               string
``pointindices``       indices into MMMx:atom atom arrays for CA/C4`   int
``populations``        populations for all *R* rotamers                (*R*,1) double
``sheet``              DSSP information on sheets                      (1,2) double
====================== =============================================== ================================ 

DSSP information (``dssp``, ``sheet``) exists only if DSSP was run or accessed through ChimeraX. Otherwise, empty output is returned. 
Furthermore, DSSP assignments refer to only the first conformer. Use ``get_conformer`` or ``cx_get_conformer`` for ensemble analysis.
Depending on the task, ``get_conformer`` may be more convenient for obtaining DSSP information even if you operate on a single chain in a single conformer. 

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
``ecoor``              argout{1} extended coordinates, argout{2}
                       corresponding atom indices
``ecoor, argout{1}``   column 1: atomic numbers; 2-4: Cartesian        (*A*,5) double
                       coordinates, 5: populations
``ecoor, argout{2}``   atom indices on entity level                    (*A*,1) int
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

-----------------------------

.. _get_conformer:
	 
Single conformer 
------------------------------------

The function ``get_conformer`` is special in that it accesses a single conformer in the entity by number.
General addressing and simultaneous access to several conformers are not supported. The function is intended for ensemble analysis.

There is not pendant ``set_conformer``. However, you can modify the returned coordinates ``coor.xyz`` and reassign them by the call 
``entity = set_coor(entity,coor.xyz,coor.indices)``.

.. code-block:: matlab

    argout = get_conformer(entity,attribute)
    [argout,exceptions] = get_conformer(entity,attribute)
    argout = get_conformer(entity,attribute,conformer)
    [argout,exceptions] = get_residue(entity,attribute,conformer)


Parameters
    *   ``entity`` - entity in MMMx:atomic format
    *   ``attribute`` - see table below (string)
    *   ``conformer`` - MMMx:atomic format conformer number, defaults to 1
Returns
    *   ``argout`` - output arguments (*M*-element cell array)
    *   ``exceptions`` - error message, if attribute is not supported  (1-element cell array)
	
**Attributes**
	
========================== =================================================== ================================
Variable                   Explanation                                         Type   
========================== =================================================== ================================
``dssp``                   DSSP information for *C* peptide chains             (1,C) array of struct
``dssp(c).chain``          chain identifier for chain *c*                      string
``dssp(c).resnum``         residue numbers                                     (1,R) int
``dssp(c).sequence``       amino acid sequence for chain *c*                   string
``dssp(c).secondary``      secondary structure assignment for chain *c*        string
``dssp(c).sheets``         sheet information for *R* residues of chain *c*     (R,2) double
``dssp(c).bp``             bridge partners for *R* residues of chain *c*       (R,2) double
``dssp(c).acc``            water exposed surface for *R* residues (Å^2)        (1,R) double
``dssp(c).NHO``            N-H-->O hydrogen bonds (two) with subfields         (1,R) array of struct
``dssp(c).NHO(r).hp``      hydrogen bonding partner for residue *r*            int
``dssp(c).NHO(r).en``      energy (kcal/mol) for residue *r*                   double
``dssp(c).OHN``            O-->H-N hydrogen bonds (two) with subfields         (1,R) array of struct
``dssp(c).OHN(r).hp``      hydrogen bonding partner for residue *r*            int
``dssp(c).OHN(r).en``      energy (kcal/mol) for residue *r*                   double
``dssp(c).tco``            cosine of angle between C=O of residues *r*, *r* +1 (1,R) double
``dssp(c).kappa``          angle between CA of residues *r*-2, *r*, *r* +2     (1,R) double
``dssp(c).alpha``          CA torsion of residues *r*-1, *r*, *r* +1, *r* +2   (1,R) double
``dssp(c).phi``            backbone dihedral phi                               (1,R) double
``dssp(c).psi``            backbone dihedral psi                               (1,R) double
``coor``                   coordinates, elements, indices                      struct
``coor.xyz``               Cartesian coordinates of all *N* atoms              (N,3) double
``coor.elements``          atomic numbers of all *N* atoms                     (N,3) double
``coor.indices``           atomic indices of *N* atoms into ``entity.xyz``     (1,N) int
``coor.all_indices``       full indices in ``entity.index_array`` format       (N,5) int16
========================== =================================================== ================================ 


Coordinates & atomic numbers (any level)
------------------------------------------

Retrieval of just the Cartesian coordinates (and optionally of the atomic numbers) is much faster with ``get_coor``.
For instance, selection of all atom coordinates in all chains of light harvesting complex LHCII (PDB 2bhw) speeds up by a factor of 17
when using ``get_coor`` instead of ``get_locations`` with attribute ``coor``. 
Unlike the object-oriented ``get`` functions, ``get_coor`` expands *all selections on different hierarchy levels* down to location level.

.. code-block:: matlab

    coor = get_coor(entity)
    coor = get_coor(entity,address)
    coor = get_coor(entity,address,heavy)
    coor = get_coor(entity,address,heavy,paradigm)
    [coor,indices] = get_coor(entity[,address[,heavy[,paradigm]]])
    [coor,indices,exceptions] = get_coor(entity[,address[,heavy[,paradigm]]])
    [coor,indices,exceptions,elements] = get_coor(entity[,address[,heavy[,paradigm]]])


Parameters
    *   ``entity`` - entity in MMMx:atomic format
    *   ``address`` - MMMx address for object selection, 'selected' or empty uses current selection
    *   ``heavy`` - flag, if true, hydrogen atoms are neglected, defaults to false
    *   ``paradigm`` - flag, if true, only the first loction/rotamer in the first conformer is returned for each atom
Returns
    *   ``coor`` - Cartesian coordinates, (*N*,3) double array for *N* selected atom locations
    *   ``indices`` - indices into entity atom tables, (*N*,1) int array
    *   ``exceptions`` - cell array of MException objects that occurred upon selection by address
    *   ``elements`` - atomic numbers, (*N*,1) int8 array

