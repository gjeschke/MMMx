Functional entity
======================

The object of interest
----------------------

Structural biology is about macromolecules, such as proteins and nucleic acids, and their complexes with other macromolecules or small ligands.

Such an object is a functional biological entity, as is assembled for fulfilling a certain task in a cell. In MMMx, we briefly call it *entity*.

Only very rarely, entities have structures, where each atom position is defined on a length scale corresponding to a typical bond length. Such high order would conflict with the dynamic nature of their function.

Usually, an entity is disordered to different extent on different levels of structural hierarchy and along the sequence of a macromolecular chain.

Description of the structure of an entity must take this into account. This requires an ensemble description.

Structural hierarchy
--------------------------------

An entity has several levels of structure

	- Primary structure
		Primary structure is fixed within an entity

			- chemical structure of each individual component
			- defines topology of the structure
			- has no explicit relation to three-dimensional structure
			- fixes some aspects of three-dimensional structure with high precision (bond lengths and bond angles)
			
	- Secondary structure
		Local conformation of a molecule 

			- often used as a synonym for the presence of secondary structure elements, such as helices and strands
			- the term is somewhat loose: three-dimensional form of *local segments* of proteins
			- there are two contributions: local backbone conformation and sidechain rotamers
			- sidechain rotamers are not usually understood as secondary structure
			
	- Tertiary structure
		Global conformation of a single molecule 

			- overall three-dimensional arrangement of a polypeptide or polynucleotide chain
			- conformation of small-molecule ligands is analogous, but not usually considered as tertiary structure
			
	- Quaternary structure
		Number and arrangement of the components of a complex 

			- while the number of components my be variable, we keep it fixed for an entity
			- different relative arrangements of subunits are usually related to different teriary structure of the subunits

----------------------------------------

This hierarchy is useful for understanding some aspects of structure, but not necessarily for representing structure of an entity or analyzing all its aspects.