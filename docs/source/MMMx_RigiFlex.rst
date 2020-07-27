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

Each rigid body has its separate coordinate array, corresponding to its *reference position and orientation*. 
In Matlab, these arrays are organized in a cell vector, which is a field `rigidbodies` of the entity variable.

Rigid bodies can be named. For this, there is a cell vector of strings `rigidbody\_names` at entity level. 

The topology variable has an additional field `rigidbody` at residue level. The value is zero for residues in flexible sections.
Otherwise, it is the number of the rigid body, which is also the index into the cell vectors `rigidbodies` and `rigidbody\_names`.

The rigid-body arrangement information stored at entity level in a field `rba`. 
For an entity with `c` conformers, `rba` is a `c \times 6` array `[\Delta x, \Delta y, \Delta z, \alpha, \beta, \gamma]` 
where `\Delta x`, `\Delta y`, and `\Delta z` are Cartesian translations (in Angstroems) and `\alpha`, `\beta`, and `\gamma` are Euler angles (`zyz` convention) in radians.

Conversion to MMMx|atomic representation
----------------------------------------

In order to avoid information loss, the RigiFlex representation is converted to the atomic representation only on two occasions

* a method that can only process the atomic representation requests the entity

* a PDB representation is required for file saving or transmission to ChimeraX

In the latter case, the MMMx|atomic representation is only used in further expansion to the PDB representation.

Expansion from `MMMx|RigiFlex` to `MMMx|atomic` explicitly creates the atom coordinates for all rigid-body arrangements, 
appends them to coordinate array `\mathbf{C}` and correspondingly amends the index array `\mathbf{I}_\mathrm{c}`. 

The fields `rigidbodies` and `rigidbody\_names` at entity level and `rigidbody` at residue level are removed only upon request. 
They allow for later back conversion, if, to a sufficiently good approximation, the domains moved as rigid bodies during 
processing by the method that requested the atomic representation.