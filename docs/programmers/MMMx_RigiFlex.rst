.. _MMMx_RigiFlex:

MMMx|RigiFlex
==========================

Basic considerations 
-----------------------------

Quite a few proteins contain *Anfinsen domains* that are ordered on the length scale of a fraction of a bond length. Usually these are small globular domains.

Representing such proteins in the `PDB` or `MMMx|atomic` format has two disadvantages

* information on rigid-body character of these domains is lost

* in large ensembles, unneccessary storage of many coordinates related by the same translation and rotation is wasteful

As MMMx creates *RigiFlex* ensemble structures, where chains can be rigid bodies joint 
by flexible linkers and rigid-body arrangement within a chain or between chains can be distributed, it makes sense to have a special representation.

Format
-------

Flexible sections are treated as in the `MMMx|atomic` format.

The coordinates for atoms in rigid bodies are stored only once, corresponding to the *reference position and orientation* of the rigid body. 
In addition, the entity variable holds index tables for all rigid bodies in a cell vector `rigidbodies`.

Rigid bodies can be named. For this, there is a cell vector of strings `rigidbody\_names` at entity level. The assignment of chains to rigid bodies is stored at entity level in a cell string `rba\_chains`.

The topology variable has an additional field `rigidbody` at chain level. It is empty for residues in flexible sections.
Otherwise, it is the number of the rigid body, which is also the index into the cell vectors `rigidbodies` and `rigidbody\_names`.

The rigid-body arrangement information stored at entity level in a field `rba`. 
For an entity with `c` conformers, `rba` is a `c \times n_{rb} \times 6` array `[\Delta x, \Delta y, \Delta z, \alpha, \beta, \gamma]` 
where `\Delta x`, `\Delta y`, and `\Delta z` are Cartesian translations (in Angstroems) and `\alpha`, `\beta`, and `\gamma` are Euler angles (`zyz` convention) in radians and `n_{rb}` is the number of rigid bodies.
Rigid-body arrangements can be assigned populations, which are stored in field `rba\_populations` at entity level.

Conversion to MMMx|atomic representation
----------------------------------------

In order to avoid information loss, the RigiFlex representation is converted to the atomic representation only on two occasions

* a method that can only process the atomic representation requests the entity

* a PDB representation is required for file saving or transmission to ChimeraX

In the latter case, the MMMx|atomic representation is only used in further expansion to the PDB representation.

It is sufficuient and much more memory efficient to sequentially request atomic representations of only individual rigid-body arrangements.
This is handled by the function ``get_rba`` which transforms atom coordinates for all rigid bodies.


