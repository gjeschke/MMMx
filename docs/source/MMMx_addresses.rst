.. _MMMx_addresses:

Selection by address
==========================

Topology and ensemble addressing
---------------------------------

Addresses are used for selecting part of an entity for display, processing, or analysis. 
Such selection can refer to objects in the molecular structure (chains or molecules in a complex, residues, atoms) or to different realizations of an object (conformers of the whole structure, sidechain rotamers, atom locations).
The first selection type is *topology addressing*  and the second selection type *ensemble addressing*.

MMMx internally uses only one index for atom locations (old style) and rotamers (new style). 
Since atom locations in crystal structure cannot generally and umabigously assigned to rotamers, addressing allows reference to either locations or rotamers in a given structure.
Mixing both concepts is discouraged, but the address style and selection functions do not prevent it.
MMMx stores occupation of water atoms from the PDB file, but does not allow to select water locations by address.

Chains
--------

Chain identifiers should be upper-case letters for up to 26 chains and lower-case letters for chains 27-52 chains. Chain identifiers longer than on letter are not compatible with PDB format, but are allowed in MMMx if required.
Chain identifiers in MMMx must start with a letter and must not contain an underscore ``_``. Chain addresses are enclosed in parentheses (*chain*). 

.. admonition:: Chain address examples

     ``(C)``  addresses chain C
	 
     ``(A,C,g)`` addresses chains A, C, and g
	 
     ``(*)`` addresses all chains
	 
   
Residues
---------

Residues can be addressed either by number or by type. Residue numbers are positive and do not neccessarily agree with the ones in the PDB file. Whether they do, is indicated in field ``original_residue_numbers`` at entity level.

Residue types are three-letter codes for amino acid residues, two-letter codes for DNA nucleotides, and single-letter codes for RNA nucleotides. Defined cofactors usually have three-letter codes. In MMMx as in ChimeraX, residue types are case-insensitive.
 
.. admonition:: Residue address examples

     ``(A)131``  addresses residue 131 in chain A
	 
     ``37-39,55`` addresses residues 37, 38, 39, and 55 in all chains
	 
     ``(C)arg,lys`` addresses all arginine and lysine residues in chain C
	 

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

Water molecules
----------------

The address ``water`` can be used to select or unselect all water molecules. Selection of individual water molecules or water atom locations by address is not supported.

Conformers
----------

Conformers are addressed by numbers in curly brackets ``{}``. By default, conformer 1 (model 1 in PDB) is selected. It is not possible to unselect all conformers. 
The selection is empty if no chain, no residue, and no atom is selected.

.. admonition:: Conformer address examples

     ``{3-5,12}``  addresses conformers 3, 4, 5, and 12 of the whole entity
	 
     ``{7}(A).CA`` addresses all ``CA`` atoms of chain A in conformer 7
	 
     ``{5}cys,met`` addresses all cysteine and methionine residues in conformer 5 of the entity

Rotamers
----------

Rotamers are addressed at residue level by their numbers after a vertical bar ``|``. By default, rotamer 1 is selected. It is not possible to unselect all rotamers if the residue, the chain, or the conformer is selected.

.. admonition:: Rotamer address examples

     ``{3-5,12}``  addresses conformers 3, 4, 5, and 12 of the whole entity
	 
     ``{7}(A).CA`` addresses all ``CA`` atoms of chain A in conformer 7
	 
     ``{5}cys,met`` addresses all cysteine and methionine residues in conformer 5 of the entity


Conformation
------------

The only information on conformation that the topology variable needs to hold are population vectors (vide supra) on entity level and on residue level (rotamers).
The length of these vectors specifies the number of entity conformers or sidechain rotamers. The topology variable must also hold the topological indices (vide supra). 
With that information, all atom coordinates can be retrieved via the index array from the coordinate array and the associated populations can be computed.

Element assignment of atoms
---------------------------

The information still missing even for basic computations concerns element assignment of the atoms, as retrieval from the atom name is not safe. 
For ease of access, this information is stored in a separate element vector that matches the coordinate index array, but can be typed `int8`. 

Conversion to PDB representation
--------------------------------

In the ``MMMx|atomic`` representation, the same atom coordinate can apply in several conformers (*models* in the PDB representation). 
Upon conversion to PDB, the coordinate array expands. The PDB writer of MMMx expands per conformer during writing to reduce memory requirements.

Rotameric states are expressed by alternate atom locations. Up to 26 (preferably) or 52 (with lower-case location identifiers) rotamers can be converted.
Not all external programs may be able to process PDB files with more than 26 locations. By default, only the 26 rotamers with highest populations are converted.
As an option, 52 rotamers can be converted.

If an atom coordinate in the ``MMMx|atomic`` representation is "not a number", this atom is ignored. 
This should happen only if the structure originated from an inconsistent PDB file.

MMMx converts to PDB representation only for two purposes:

* saving structure in a PDB files

* transmitting structure to ChimeraX for visualization

Conversion from PDB representation
----------------------------------

MMMx does not make an effort to preserve atom numbers and only a limited effort to preserve residue numbers of the original PDB entry. 
Residue numbers are preserved in the about `96.5% structures that do not use "insertion codes"`__ and only if all are positive numbers.  
The entity has a field ``original_element_numbers`` that indicates whether residue numbers were preserved.

.. __: http://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v40.dic/Items/_atom_site.pdbx_PDB_ins_code.html

The number of chain/molecule conformers equals the number of PDB ``models`` for all chains and molecules of the entity. Uniform populations are assumed.

The number of rotamers of a residue is as large as the maximum number of alternate locations among the atoms of this residue. 
Rotamer populations are mean populations over all atoms which have this number of alternate locations.

In case of topological inconsistency between *models*, topology is determined by the first model encountered in the PDB file (regardless of its model number).
Only atom coordinates are read for further models. Surplus atom coordinates are ignored. Missing atom coordinates are assigned "not a number".