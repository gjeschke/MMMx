ChimeraX
======================


MMMx communicates with ChimeraX for visualization. MMMx can use ChimeraX for PDB and mmCIF file loading and saving, for retrieving information on an entity, for entity building, and entity modification.

The graphical user interfaces (GUIs) of ChimeraX and MMMx can be used simultaneously. 
The master copy of an entity is hold by MMMx.   

	- Communication protocol
		MMMx uses REST (representational state transfer) remote control both for setting states of ChimeraX and for polling ChimeraX for information:

			- ChimeraX needs to call the command ``remotecontrol rest start port 51051`` upon startup
			- communication depends on log response of ChimeraX commands
			- MMMx claims the name of entities (models) that it loads into ChimeraX

	- Using MMMx without ChimeraX
		MMMx can be used standalone, including its own GUI, but the following restrictions apply: 

			- no 3D structure visualization
			- no mmCIF support (will be changed later)
			- not all selection modes of ChimeraX are supported
			- for methods that require secondary structure analysis, DSSP must be installed and registered

-----------------------

Remarks
-----------------------

* The startup commands can be set in ChimeraX via the menu item ``Favorites\Settings...`` after navigating to tab ``Startup``.

* The log response of ChimeraX is very well structured, but we are not guarded against future format changes

* MMMx models in ChimeraX have unique names MMMx_*mmmid*, where *mmmid* is the internal identifier in MMMx
  This guards against accidental replacement of an MMMx model by another model by the user. 
  The other model could then have the same number, but would not have the MMMx name.

* If ChimeraX is running and connected, MMMx can access DSSP secondary structure analysis through ChimeraX

