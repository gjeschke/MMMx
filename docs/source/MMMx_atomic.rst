.. _MMMx_atomic:

MMMx|atomic
==========================

Basic considerations
-----------------------------

From a thermodynamic point of view, different conformers have different probability to occur, as they have different free energies. They should thus be assigned populations.

If modelling is based on restraints in the form of probability density distributions, such as distance distribution restraints, probabilities of conformers are accessible in modelling.

The reason for alternate atom locations in protein crystal structures is, in most cases, significant population of more than one sidechain rotamer. 
This should be assigned at residue level, not atom level.

From both a conceptual and data processing point of view, it makes sense to separate topology (fixed and the same for all conformers) from coordinates (variable).

Coordinates are heterogeneously different, meaning that for some atoms they do differ between conformers, while for others they do not. 
Especially for large ensembles, it makes sense to store only unique coordinates.

.. admonition:: Coordinate indices

     An atom coordinate is assigned by five indices. Three are *topological* and designate
	 
	 1) the molecule, which can be a macromolecular chain
	 
	 2) the residue in a chain
	 
	 3) the atom in a residue
	 
	 The remaining two indices are *conformational* and designate
	 
	 4) the conformer
	 
	 5) the rotamer of a sidegroup (structures with specified rotamers) or one of the atom locations existing for a residue (crystal structures with locations)
   
The Cartesian coordinate array :math:`\mathbf{C}` of an entity is thus supplemented by a coordinate index array :math:`\mathbf{I}_\mathrm{c}`. 
In Matlab, `logical indexing`__, together with the index array, provides a very intuitive and efficient way for extracting required subsets of coordinates and to reassign them after some processing. 
The index array is typed ``uint16``, so that each hierarchy level can hold up to 65535 items.

.. __: https://blogs.mathworks.com/loren/2013/02/20/logical-indexing-multiple-conditions/

The coordinate array is further supplemented by a vector of atomic numbers (number of the element in the periodic table) and by a vector of occupancies. Both vectors are single-byte unsigned integers. Occupancy is multiplied by 100 and rounded to an integer for storage.

Topology
----------

Topology describes chemical structure and is immutable for an entity. In other words, any chemical reaction, such as posttranslational modification, is a transfer from one entity to another entity.

Topology is represented in terms of the number of objects (chains, molecules, residues, atoms), their identifiers, and their hierarchical relation. 
Hierarchical relation means that an atom belongs to a residue or molecule and a residue belongs to a chain. The term *residue* includes amino acid residues, nucleotide residues, and, if someone finds a way to characterize carbohydrates, a sugar residue.

In Matlab, hierarchy is easily represented by a multi-level structure. By combining this with `dynamic fieldnames`__, storage and access become very convenient and efficient.

.. __: https://ch.mathworks.com/help/matlab/matlab_prog/generate-field-names-from-variables.html

Further attributes of objects can be stored by extending levels of the topology variable with additional fields.

.. admonition:: Field name conventions

     All object names start with capital letters. They are:
	 
	 - **Chain identifiers**. Small letters from PDB files with more than 26 chains are represented by a capital letter with an appended underscore. 
	   Entities with more than 26 chains can be exchanged with other software only in ``mmCIF`` format.
	   
	 - **Residue identifiers**. Capital ``R`` followed by the residue number.
	 
	 - **Atom identifiers**. The atom name in PDB or mmCIF syntax. Primes are substituted by an underscore.
	   If a name starting with a lowercase character is encountered, the first letter is capitalized and an underscore is appended to the atom name.
	   
	 The fieldname ``selected`` is reserved at all levels. At entity level it is a vector of selected conformers (default: 1), at other levels a boolean flag.
	 
	 The fieldname ``index`` is reserved at all levels below entity. This field holds the topology index per object and level.
	 
	 The fieldnames ``water`` and ``water_selected`` at entity level are reserved for water molecules. They hold an index vector of all water atoms and a flag whether they are selected.
	 
	 The fieldnames ``name``, ``xyz``, ``index_array``, ``elements``, ``occupancies``, and ``original_residue_numbers`` are reserved at entity level

	 The fieldname ``populations`` is reserved at entity level and at residue level. At residue level, the length of the population vecor specifies the number of rotamers.
	 
	 The fieldnames ``name``, ``selected_rotamers``, ``locations``, ```dssp``, and ``sheet`` are reserved at residue level
	 
	 The fieldnames ``element``, ``selected_locations``, ``charge``, and ``bfactor`` are reserved at atom level
	 
	 Methods can define further lower-case fieldnames for object attributes. These attributes are internal to the method. 
	 Passing back such attributes to MMMx and saving them in ``mmCIF`` output files may be enabled in a future release. 


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

Upon loading a PDB file, MMMx does not make an effort to preserve atom numbers and only a limited effort to preserve residue numbers of the original PDB entry. 
Residue numbers are preserved in the about `96.5% structures that do not use "insertion codes"`__ and only if they all are positive numbers and if, within the same chain, they appear in ascending order in the PDB file.  
The entity has a field ``original_residue_numbers`` that indicates whether residue numbers were preserved.

.. __: http://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v40.dic/Items/_atom_site.pdbx_PDB_ins_code.html

The number of chain/molecule conformers equals the number of PDB ``models`` for all chains and molecules of the entity. Uniform populations are assumed.

The number of rotamers of a residue is as large as the maximum number of alternate locations among the atoms of this residue. 
Rotamer populations are mean populations over all atoms which have this number of alternate locations.

In case of topological inconsistency between *models*, topology is determined by the first model encountered in the PDB file (regardless of its model number).
Only atom coordinates are read for further models. Surplus atom coordinates are ignored. Missing atom coordinates are assigned "not a number".