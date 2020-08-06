.. _distance_distributions:

Distance distributions
==========================

Concept
---------------------------------

Methods in MMMx should not directly change the ``entity`` variable, but rather use the ``set_*object*`` functions, where *object* is a conformer, a chain, a residue, a rotamer, an atom, or an atom location. 
The ``entity`` is an MMMx ensemble structure representation in :ref:`MMMx|atomic<MMMx_atomic>` or :ref:`MMMx|RigiFlex<MMMx_RigiFlex>` format.

``set`` functions modify attributes of objects, whereas ``get`` functions retrieve them. They are specific to an object hierarchy level in order to avoid unintended behaviour.
The objects are selected by an :ref:`MMMx address<MMMx_addresses>`. The address can also be ``selected`` for accessing all currently selected objects on this hierarchy level. 
If an address is provided, selection in the entity changes to the one specified by this address.

It is good practice to use ``get_"object"`` to retrieve current attributes, modify them, and reassign the modified attributes by ``set_object``.
Between the calls of ``get_"object"`` and ``set_"object"``, either the selection in the entity should not be changed (calls without ``address``)
or the same address must be used in both calls. 

The special function ``set_coor`` can be used to assign modified coordinates. 
It can be applied to a mixed selection on several hierarchy levels and is faster than coordinate access by the ``set_"object"`` functions.

For transforming coordinates of all atoms of a single conformer, the fastest way is to retrieve coordinates and atom indices with :ref:`get_conformer<get_conformer>` 
and to use ``entity = set_coor(entity,coor.xyz,coor.indices)`` for reassignment of modified coordinates. For changing conformer populations, directly access ``entity.populations``. 

The attribute ``info`` cannot generally be modified, as this could compromise integrity of addressing or indexing the entity. 
There are exceptions for fields of info on some hierarchy levels.

Generic syntax
--------------

.. code-block:: matlab

    [entity,exceptions] = set_"object"(entity,attribute,argin,address)
	 
where the ``attribute`` strings are specific to the "object" hierarchy level 
("residue", "atom", "location", see below) and argin is a cell vector, 
whose length is (most of the time) the number of objects of this hierarchy level that were selected by ``address``.
If the address is ``selected`` or empty, the ``set`` functions operate on objects currently selected in the entity.

There is no ``set`` function for chains. Chain tags can be changed only upon saving the entity, for instance with ``put_pdb``.

Error messages or warnings are reported as MException objects in cell array ``exceptions``. 

 
Residue
---------

Use this function if you want to operate exclusively on residue level.
Selections above or below residue level are ignored.

.. code-block:: matlab

    entity = set_residue(entity,attribute,argin)
    [entity,exceptions] = set_residue(entity,attribute,argin)
    entity = set_residue(entity,attribute,argin,address)
    [entity,exceptions] = set_residue(entity,attribute,argin,address)

Parameters
    *   ``entity`` - entity in MMMx:atomic format
    *   ``attribute`` - see table below (string)
    *   ``address`` - MMMx address for object selection, 'selected' or empty uses current selection
Returns
    *   ``entity`` - modified entity in MMMx:atomic format
    *   ``exceptions`` - error message, if attribute is not supported  (1-element cell array)
	
**Attributes**
	
====================== =============================================== ================================
Variable               Explanation                                     Type   
====================== =============================================== ================================
``dssp``               DSSP secondary structure assignment             char
``number``             residue number                                  int
``populations``        populations for all *R* rotamers                (*R*,1) double
``sheet``              DSSP information on sheets                      (1,2) double
``tlc``                three-letter code/PDB residue tag               string
====================== =============================================== ================================ 

``set_residue`` raises an exception and returns an empty entity if the requested residue number already exists in that chain.

Three-letter codes are capitalized, any characters that are not alphanumeric are removed, and the string is trimmed and truncated to at most three characters. 
If the request does not contain alphanumeric characters, an exception is raised and the returned entity is empty. 

The format of ``dssp`` and ``sheet`` information is not checked.

-----------------------------

Atoms
---------

Use this function if you want to operate on atoms of all rotamers or on all atom locations.
Selections above atom level are ignored.

Use Matlab built-in function ``cell2mat`` for reforming output for B factor, charge, atomic number, and population into vectors. 
Note that MMMx supports only one B factor per atom, not distinct B factors for locations.

.. code-block:: matlab

    entity = set_atom(entity,attribute,argin)
    [entity,exceptions] = set_atom(entity,attribute,argin)
    entity = set_atom(entity,attribute,argin,address)
    [entity,exceptions] = set_atom(entity,attribute,argin,address)

Parameters
    *   ``entity`` - entity in MMMx:atomic format
    *   ``attribute`` - see table below (string)
    *   ``argin`` - input arguments (*M*-element cell array for *M* selected locations)
    *   ``address`` - MMMx address for object selection, 'selected' or empty uses current selection
Returns
    *   ``entity`` - modified entity in MMMx:atomic format
    *   ``exceptions`` - error message, if attribute is not supported  (1-element cell array)
	
**Attributes**
	
====================== =============================================== ================================
Variable               Explanation                                     Type   
====================== =============================================== ================================
``bfactor``            crystallographic B factor                       double
``charge``             atom charge                                     int
``coor``               Cartesian coordinate array for *all* locations  (*N*,3) double
``element``            atomic number                                   int8      
``name``               atom name                                       string, maximum 4 characters  
``population``         rotamer population or atom occupancy            double
``xyz``                Cartesian coordinates per location              (1,3) double
====================== =============================================== ================================ 

-----------------------------

Modifying atom names is generally discouraged, but may be useful for (paramagnetic) substitution of ions.
The atom name is capitalized, primes are substituted by underscores, and it is truncated to 4 characters.
   
Locations
---------

Use this function if you want to operate on selected rotamers or atom locations.
If the selection is on atom level and no rotamers are selected, only the first location or rotamer is referred to.
Selections above atom level are ignored.

.. code-block:: matlab

    entity = set_location(entity,attribute,argin)
    [entity,exceptions] = set_location(entity,attribute,argin)
    [entity,exceptions] = set_location(entity,attribute,argin,address)


Parameters
    *   ``entity`` - entity in MMMx:atomic format
    *   ``attribute`` - see table below (string)
    *   ``argin`` - input arguments (*M*-element cell array for *M* selected locations)
    *   ``address`` - MMMx address for object selection, 'selected' or empty uses current selection
Returns
    *   ``entity`` - modified entity in MMMx:atomic format
    *   ``exceptions`` - error message, if attribute is not supported  (1-element cell array)
	
**Attributes**
	
====================== =============================================== ================================
Variable               Explanation                                     Type   
====================== =============================================== ================================
``element``            atomic number                                   int8        
``population``         rotamer population or atom occupancy            double
``xyz``                Cartesian coordinates per location              (1,3) double
====================== =============================================== ================================ 

-----------------------------
	 

Coordinates (any level)
------------------------------------------

For modification of only Cartesian coordinates of a set of objects, it is faster to use ``get_coor`` and ``set_coor``.

If no atom indices are provided, ``set_coor`` expands *all selections on different hierarchy levels* down to location level.
The function is much faster when the atom indices, originally retrieved by ``get_coor``, are provided instead of an ``address``.
However, in order to enable better code readability, a call with ``address`` or for the current selection is allowed.

.. code-block:: matlab

    [entity,exceptions] = set_coor(entity,coor)
    [entity,exceptions] = set_coor(entity,coor,indices)
    [entity,exceptions] = set_coor(entity,coor,address)

    
Parameters
    *   ``entity`` - entity in MMMx:atomic format
    *   ``coor`` - Cartesian coordinates, (*N*,3) double array for *N* selected or indexed atom locations
    *   ``indices`` - indices into entity atom tables, (*N*,1) int array
    *   ``address`` - MMMx address for object selection, 'selected' or empty uses current selection, string
Returns
    *   ``entity`` - modified entity in MMMx:atomic format
    *   ``exceptions`` - cell array of MException objects that occurred upon selection by address

