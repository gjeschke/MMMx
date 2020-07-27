Data
======================

Entity
------

An *entity* holds the structure information on a protein, a nucleic acid molecule, or complexes thereof, possibly including small-molecule ligands and water molecules.  

The entity consists of several *conformers* for an ensemble structure. Conformers have the same primary structure (amino acid or nucleotide sequence), but can differ in secondary, tertiary, and quaternary structure.

Each atom has a unique address and, in general, Cartesian coordinates (in Angstroem). Cartesian coordinates can be unspecified ("not a number").

.. admonition:: Addresses

     The atom address for specifying the CA atom of residue 131 in conformer 1 of chain A in structure 2LZM is **``[2LZM](A){1}131.CA``** in MMM. In ChimeraX, it is **``#1/A:131@CA``**. The user must know that structure 2LZM is model #1.

Methods can operate on the coordinates, which are hold in one common array for each conformer of a chain. Methods can also change residues and thus sequence in a chain (mutation). For any changes beyond that, a new entity should be derived. 

.. admonition:: Matlab realization

     In Matlab, an entity is described by data type ``struct`` with dynamic field names. Field names that begin with a capital letter describe primary structure (chain identifiers, residue numbers, conformer numbers, atom identifiers). Field names that begin with a lower-case letter provide additional information, such as coordinates, residue names, the element, etc.

The MMMx entity file format is based on this extensible Matlab structure variable. Extension is by *attributes*. Each level (entity, chain, residue, conformer, atom) can have additional attributes, whose lower-case names are dynamical field names in Matlab.

Rotamer library
--------------------------------

A rotamer library contains sidechain conformers for a spin label or chromophore label, together with their populations in the absence of interactions with a macromolecule.

Further *attributes* specify the attachment frame (by three atoms), the label position and possibly the label orientation. They may also provide information on how the library was generated.

Representation in Matlab and as MMMx file format is analogous to an *entity*.


Distance distribution restraints
--------------------------------

Distance distribution restraints are specified at least by a distance axis (in Angstroem) and probabilities that a distance falls into a given bin on the diatance axis. The sum of all probabilities is unity.

Further, optional attributes are lower and upper probability bounds, a handle to a function that generates a parametric distance distribution, and parameter values and uncertainties.

