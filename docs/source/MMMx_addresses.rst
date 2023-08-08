.. _MMMx_addresses:

Selection by address
==========================

Topology and ensemble addressing
---------------------------------

Addresses are used for selecting part of an entity (protein, nucleic acid, or complex). 
Such selection can refer to objects in the molecular structure (chains or molecules in a complex, residues, atoms)
or to different realizations of an object (conformers of the whole structure, atom locations).
The first selection type is *topology addressing*  and the second selection type *ensemble addressing*.

MMMx internally uses only one index for atom locations (old style) and rotamers (new style). 
MMMx stores occupation of water atoms from the PDB file, but does not allow to select water locations by address.

Chains
--------

Chain identifiers should be upper-case letters for up to 26 chains and lower-case letters for chains 27-52. 
Chain identifiers longer than on letter are not compatible with PDB format and not allowed in MMMx.
Chain addresses are enclosed in parentheses (*chain*). 

.. admonition:: Chain address examples

     ``(C)``  addresses chain C
	 
     ``(A,C,g)`` addresses chains A, C, and g
	 
     ``(*)`` addresses all chains
	 
   
Residues
---------

Residues can be addressed either by number or by type. Residue numbers are positive and do not neccessarily agree with the ones in the PDB file. 
Whether they do, is internally indicated in field ``original_residue_numbers`` at entity level. The first residue
in MMMx has number 1 if the PDB file contains negative residue numbers or the residue number zero. If the PDB file starts 
with a residue number larger than 1, residue numbers do agree. The only exception to these rules are PDB structures that use the infamous 
insertion code. MMMx counts up, so that all residues from the inserted residue on have numbers increased by one.

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

The address ``water`` can be used to select or unselect all water molecules. 
Selection of individual water molecules or water atom locations by address is not supported.

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

Rotamer addressing is not currently used in modeling modules.

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
	 

Changes compared to MMM
-----------------------

The MMMx address format was designed to be as close as possible to the MMM address format, yet allowing for full access to the ``MMMx:atomic`` representation of ensemble structure.
This entailed the following changes:

* the structure identifier in square brackets is no longer required, since MMMx methods work on an entity

* rotamer addressing was newly introduced

* the wildcard is now the asterisk ``*`` rather than the colon ``:``

* preferably, the conformer is now addressed first, whereas MMM addressed it as "chain model" after the chain identifier; both address sequences still work

Correspondence with ChimeraX
----------------------------

MMMx has a function that can translate a subset of basic ChimeraX target specifications into MMMx addresses. The following functionality of ChimeraX target specifications is **not** translated:

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

MMMx can also project its own addresses onto ChimeraX target specifications, as far as the ChimeraX target specification supports the type of addressing. 
This excludes rotamer and atom location addressing.
