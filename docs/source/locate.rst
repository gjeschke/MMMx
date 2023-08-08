.. _locate:

Locate
==========================

This module computes a probability density that locates a spectroscopic label with respect to a resolved part of a structure.

Locate can have an input argument:

.. code-block:: matlab

    !locate [coverage]

Arguments
    *   ``coverage`` - probability covered by an isosurface, defaults to 0.5  
Remarks
    *   the density level corresponding to this isosurface is written to the logfile

``cube``
---------------------------------

Specifies edge length of the density cube 

.. code-block:: matlab

    cube size

Arguments
    *   ``size`` - edge length of the density cube in Angstroem, cube volume is `size \times size \times size`, default is 75
Remarks
    *   the most probable location defines the cube center
    *   the largest relative probability density at a cube border is reported in the logfile
    *   if there is high relative probability density with the default values, the location is very poorly restrained

``getAlphaFold``
---------------------------------

Input of a template by reading an AlphaFold prediction from the database. 

.. code-block:: matlab

    getAlphaFold UniProtID

Arguments
    *   ``UniProtID`` - UniProt identifier of the AlphaFold prediction
Remarks
    *   note that not all proteins in UniProt have an AlphaFold prediction in the database
	
``getpdb``
---------------------------------

Input of a template by reading a single local PDB file or retrieving a PDB structure from a server. 

.. code-block:: matlab

    getpdb file

Arguments
    *   ``file`` - file name or PDB identifier, must have extension '.pdb' if it specifies a file
Remarks
    *   if the template has several conformers, the first one is used
	
``grid``
---------------------------------

Specifies number of grid points along one dimension of the density cube 

.. code-block:: matlab

    grid size

Arguments
    *   ``size`` - number of grid points along one dimension, grid size is `size \times size \times size`, default is 176
Remarks
    *   memory and computation time scale with `size^3`
	
``reference``
---------------------------------

Definition of distance distribution restraints to reference points. This is a block key with `n` lines for `n` reference points. 

.. code-block:: matlab

    reference label
       address 'rmean' 'rstd' [@'fname']
       ...
    .reference

Arguments
    *   ``label`` - label type, e.g. `mtsl`
    *   ``address`` - address of the reference point, e.g., `(A)75`
    *   ``rmean`` mean distance in Angstroem, e.g. `32.7`
    *   ``rstd`` standard deviation in Angstroem, e.g. `3.9`
    *   ``fname`` optional file name of the distance distribution 
Remarks
    *   use separate 'reference' blocks for different label types
    *   the file name is optional, full distributions can be used
    *   if a full distribution is provided, ``rmean`` and ``rstd`` can be skipped

``save``
---------------------------------

Specifies basis name for saving output conformers 

.. code-block:: matlab

    save file

Arguments
    *   ``file`` - file name for output density cube file
Remarks
    *   extension '.mat' is appended if there is none, use extension '.mat' for MMMx density files or '.mrc' for MRC files
	
