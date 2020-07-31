.. _object_access:

Object access
==========================

Concept
---------------------------------

Object access functions return or set object attributes, where an object is a conformer, a chain, a residue, a rotamer, an atom, or an atom location. 
They operate on an ``entity``, which is an MMMx ensemble structure representation in :ref:`MMMx|atomic<MMMx_atomic>` or :ref:`MMMx|RigiFlex<MMMx_RigiFlex>` format.
Wherever possible, methods should operate on an entity only through object access functions.

``get`` functions retrieve attributes, whereas ``set`` functions set or change attributes. They are specific to an object hierarchy level in order to avoid unintended behaviour.
The objects are selected by an :ref:`MMMx address<MMMx_addresses>`. The address can also be ``selected`` for accessing all currently selected objects on this hierarchy level. 

Generic syntax
--------------

.. code-block:: matlab

    [argout,exceptions] = get_"object"(entity,address,attribute)
	 
where possible ``attribute`` strings are specific to the "object" hierarchy level ("conformer", "chain", "residue", "rotamer", "atom", "location", see below) and argout is a cell vector, whose length is the number of objects of this hierarchy level that were selected by ``address``.
Error messages or warnings are reported as MException objects in cell array ``exceptions``.  
   
Residues
---------

Residues can be addressed either by number or by type. Residue numbers are positive and do not neccessarily agree with the ones in the PDB file. Whether they do, is indicated in field ``original_residue_numbers`` at entity level.

Residue types are three-letter codes for amino acid residues, two-letter codes for DNA nucleotides, and single-letter codes for RNA nucleotides. Defined cofactors usually have three-letter codes. In MMMx as in ChimeraX, residue types are case-insensitive.
 
.. admonition:: Residue address examples

     ``(A)131``  addresses residue 131 in chain A
	 
     ``37-39,55`` addresses residues 37, 38, 39, and 55 in all chains
	 
     ``(C)arg,lys`` addresses all arginine and lysine residues in chain C
	 
     ``(B)*`` addresses all residues in chain B
	 

Atoms
------

Atoms can be addressed by type, which is usually a string with up to four characters that must start with a letter. 
If an atom name in the original PDB file started with a digit, it is preceded by ``A_`` in MMMx.
The atom address ``backbone`` selects the ``N``, ``CA``, ``C``, and ``O`` atoms of a peptide residue.
A list of atom addresses starts with a dot ``.``.

.. admonition:: Atom address examples

     ``(B)22.CB``  addresses the ``CB`` atom in residue 22 of chain B
	 
     ``.CG1,CG2`` addresses all ``CG1`` and ``CG2`` atoms in the entity
	 
     ``(A)128-135.backbone,CB`` addresses the backbone atoms and the ``CB`` atom in the residue range from 128 to 135 in chain A
	 
     ``131.*`` addresses all atoms in residues 131 of all chains

Water molecules
----------------

The address ``water`` can be used to select or unselect all water molecules. Selection of individual water molecules or water atom locations by address is not supported.

Conformers
----------

Conformers are addressed by numbers in curly brackets ``{}``. By default, conformer 1 (model 1 in PDB) is selected. It is not possible to unselect all conformers. 
The selection is empty if no chain, no residue, and no atom is selected. For the conformer range, keywords ``start`` (first conformer) and ``end`` (last conformer) are supported.

.. admonition:: Conformer address examples

     ``{3-5,12}``  addresses conformers 3, 4, 5, and 12 of the whole entity
	 
     ``{4-end}`` addresses all conformers starting at number 4
	 
     ``{7}(A).CA`` addresses all ``CA`` atoms of chain A in conformer 7
	 
     ``{5}cys,met`` addresses all cysteine and methionine residues in conformer 5 of the entity
	 
     ``{*}`` selects all conformers of the entity

Rotamers
----------

Rotamers are addressed at residue level by their numbers after a vertical bar ``|``. 
By default, rotamer 1 is selected. It is not possible to unselect all rotamers if the residue, the chain, or the conformer is selected.
Rotamer selection overrules location selection for atoms. In a rotameric structure, atom locations correspond to distinct rotamers.

.. admonition:: Rotamer address examples

     ``(A)131|1-3``  addresses rotamers 1, 2, 3 of residue 131 in chain (A)
	 
     ``55-57|2.CA`` addresses the second *location* (where no location tag is often the first location) for the ``CA`` atoms of residues 55, 56, and 57 in all chains  
	 
     ``gln|2`` addresses the second rotamer for all glutamine residues in the entity
	 
     ``glu|*`` addresses all rotamers for all glutamate residues in the entity

Locations
----------

Locations are addressed at atom level by their tags after a colon ``:``.  
A location tag is preferably a single upper-case letter for up to 26 locations and a single lower-case letter for locations 27-52. More than 52 locations cannot be stored in PDB files.
By default, the first location (location tag is usually a space) is selected. It is not possible to unselect all locations if the atom is selected.
If you want to address a location by its tag, you must not address rotamers in the same address.

.. admonition:: Location address examples

     ``.OE1:B``  addresses location B of all OE1 atoms in the entity, if no such location exists for an OE1 atom, the first location is selected 
	 
     ``(A)glu.OE1:A`` addresses location A of OE1 atoms of all glutamate residues in chain A
	 
     ``glu.OE1:*`` addresses all locations of OE1 atoms in all glutamate residues in chain A 
	 
Programmatic access
-------------------

.. code-block:: matlab

    [entity,exceptions] = select(entity,address,overwrite,unselect)
    
selects objects by ``address`` in an ``entity``. If ``overwrite`` (default: false) is true, a pre-existing selection is deleted. 
If ``unselected`` (default: false) is true, the addressed objects are unselected rather than selected. An unselect request overrules a simultaneous overwrite request.

.. code-block:: matlab

    [atom_indices,complete] = get_selection(entity)
    
retrives sorted atom indices into the ``.xyz```, ``.elements``, and ``.occupancies``, and ``.index_array`` fields of ``entity``. 
If requested, ``complete`` returns full indices (chain, residue, atom, conformer, rotamer/location), which are not sorted.

.. code-block:: matlab
     
    [addresses,exceptions] = cx_to_mmmx_address(spec)	 

translates the ChimeraX target specification ``spec`` into MMMx addresses (cell array, one entry for each structure identifier in ``spec``). See below.
Error messages or warnings are reported as MException objects in cell array ``exceptions``.

.. code-block:: matlab
     
    [spec,exceptions] = cx_from_mmmx_address(address,id)

Projects an MMMx ``address`` to a ChimeraX target specification	``spec`` (see below), for which an optional structure identifier ``ìd`` in ChimeraX can be supplied.
Error messages or warnings are reported as MException objects in cell array ``exceptions``.


Changes compared to MMM
-----------------------

The MMMx address format was designed to be as close as possible to the MMM address format, but allowing for full access to the ``MMMx:atomic`` representation of ensemble structure.
This entailed the following changes:

* the structure identifier in square brackets is no longer required, since MMMx methods work on an entity

* rotamer addressing was newly introduced

* the wildcard is now the asterisk ``*`` rather than the colon ``:``

* preferably, the conformer is now addressed first, whereas MMM addressed it as "chain model" after the chain identifier; both address sequences still work

Correspondence with ChimeraX
----------------------------

MMMx can translate a subset of basic ChimeraX target specifications into MMM addresses. The following functionality of ChimeraX target specifications is **not** translated:

* usage of ``start`` and ``end`` in ranges, except for conformers
* structure hierarchy with a depth of more than two (only structure and conformers allowed)
* ranges for structure identifiers or chain identifiers
* implicit operations (returning higher-level address parts)
* use of the wild card ``*`` for part of a name
* use of the wild card ``?`` for single characters
* built-in classifications
* user-defned attributes
* zones
* combinations

Such selections can be made in ChimeraX, also via the ChimeraX interface of MMMx, and can then be imported into MMMx.

MMMx can project its own addresses onto ChimeraX target specifications, as far as the ChimeraX target specification supports the type of addressing. 
This excludes rotamer and atom location addressing.
