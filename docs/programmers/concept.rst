Concept
======================


MMMx relies on well designed and well maintained software for standard tasks.

It is intended for modelling structure of proteins and their complexes and for supporting label-based spectroscopic techniques, such as EPR, FRET and NMR techniques that use spin labels.

	- Imports
		Imports are functionalities that MMMx sources from external apps or packages:

			- Protein and nucleic acid structure loading, saving, manipulation and visualization by *ChimeraX*
			- Small-angle scattering curve simulation and fitting by the *ATSAS* package
			- Amino acid side chain generation or optimization by *SCWRL4*
			- Sequence alignment by *MUSCLE*
			- Structure refinement by *Yasara*

	- Exports
		Exports are functionalities that MMMx provides by building or processing structure data: 

			- Spin label and chromophore modelling via rotamer libraries
			- Simulation of distance distributions between labels
			- Localization of labels at structurally unresolved sites with espect to a structure
			- Conformational transitions by distance-restraint based deformation of elastic network models
			- Ensemble modelling of flexible peptide and RNA chains or chain segments
			- RigiFlex ensemble models based on distributed rigid-body arrangement and flexible segments


-----------------------

Open source concept
---------------------------------------------

All MMMx functions and definitions are open source. However, MMMx may rely for some of its functionality on external apps that are not open source or not even free (Yasara for structure refinement is an example).

Separation of user interface from computation
---------------------------------------------

The graphical user interface (GUI) of MMMx coordinates the work of MMMx methods and functionalities of external apps and visualizes numerical data. 

MMMx methods can be called from scripts without running the GUI and without reliance on ChimeraX. For that, MMMx has basic PDB file read and write facility.


Visibility to users and programmers
------------------------------------

Although MMMx is open source, users and programmers need to be aware that implementation details can change without notice.

Only part of MMMx functionality can reliably be used in scripts, contributed methods, or by external requests:

	- Data
		The following types of data have a stable, though extensible definition:

			- *entity*: internal representation of macromolecular (ensemble) structure
			- *rotlib*: rotamer library
			- *ddr*: distance distribution restraint

	- Functionality
		The following elements of MMMx functionality will be stable, though extensible: 

			- *Commands*: allow for control and scripting in the MMMx GUI
			- *Methods*: build, modify or analyze an entity
			- *Restraint file*: specify experimental restraints and methods that operate on them



 

