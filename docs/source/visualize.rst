.. _visualize:

Visualize
==========================

This module translates MMMx visualization scripts to the script language of MMM and can directly execute them in MMM. 
For direct execution, MMM must be open in the same Matlab instance. Isosurface visualization of (pseudo-)electron densities does not require MMM. 
The module is intended for convenient and consistent visualization of ensembles generated by MMMx.
Population of conformers is encoded either by transparency or by coil radius.

``addpdb``
---------------------------------

Input of an ensemble via PDB files matching a pattern. This can also be used for visualizing a single conformer 

.. code-block:: matlab

    addpdb files

Arguments
    *   ``files`` - file name that can include wildcards, e.g. ``hnRNPA1_sorted_m*.pdb``
Remarks
    *   all populations are equal
    *   each PDB file should hold only a single conformer or, alternatively, there should be only one PDB file with all conformers
    *   there should be only one ``addpdb`` line, or, alternatively, a ``getens``, an ``import``, or a ``getAlphaFold`` line 
	
``color``
---------------------------------

Set color of graphics elements. All conformers in the ensemble have the same coloring; only transparency or coil radius differs.

.. code-block:: matlab

    color address rgb

Arguments
    *   ``address`` - address of objects of a conformer, e.g. `(A)58-154` for residues 54-154 of chain A
    *   ``rgb`` - rgb specifier, either three numbers between 0 and 1 for red, green, blue or SVG color name
Remarks
    *   ``color (A)455-531 0.0 1.0 1.0`` will set a cyan color (0% red, 100% green, 100% blue)
    *   ``color (A)455-531 forestgreen`` will set SVG color `forestgreen`
    *   see `SVG color table <https://www.december.com/html/spec/colorsvg.html>`_ for available colors

``colorscheme``
---------------------------------

Set color scheme of graphics elements. All conformers in the ensemble have the same coloring; only transparency or coil radius differs.

.. code-block:: matlab

    colorscheme address scheme

Arguments
    *   ``address`` - address of objects of a conformer, e.g. `(A)58-154` for residues 54-154 of chain A
    *   ``scheme`` - one of the color schemes available in MMM, see Remarks
Remarks
    *   the following remarks explain the available schemes
    *   `secondary` - standard metal colors to secondary structure elements in ribbon plots (helices: copper, sheets: steelblue, loops: gold)
    *   `sequence` - rainbow colors from blue to red within a chain starting from the N terminus to the C terminus 
    *   `charge` - colors from blue to red according to charge. Dark blue: +2, blue +1, grey 0, red: -1, dark red: -2
    *   `hydropathy` - colors from blue to red according to hydropathy, with the blue end corresponding to the most hydrophilic and the red end to the most hydrophobic residues
    *   `helix` _ `propensity` - rainbow colors from blue to red according to helix propensity, with the blue end corresponding to residues that tend to form an alpha helix and the red end to helix-breaking residues (Pro). Color grade is proportional to the square root of helix propensity
    *   `sequence` can be used with two additional arguments, address of the first residue and number of residues for the color grade
    *   if `sequence` is used with additional arguments, all residues before and after the indicated segment are colored blue and red, respectively

``density``
---------------------------------

Add a density surface from an MMMx density file generated by module Locate or EnsembleAnalysis.

.. code-block:: matlab

    density file [level [opacity [rgb]]]

Arguments
    *   ``file`` - name of an MMMx density file, generate such files with the ``density`` keyword in EnsembleAnalysis
    *   ``level`` - fraction of total density enclosed by the isosurface, defaults to 0.5
    *   ``opacity`` - opacity of the isosurface, 1 is completely opaque, 0 is invisible, defaults to 0.5
    *   ``rgb`` - rgb specifier, either three numbers between 0 and 1 for red, green, blue or SVG color name
Remarks
    *   see `SVG color table <https://www.december.com/html/spec/colorsvg.html>`_ for available colors
    *   ``rgb`` defaults to '0.75 0 0', which is a darkish red

``execute``
---------------------------------

Requests direct execution of the visualization script in MMM.   

.. code-block:: matlab

    execute

Remarks
    *   MMM must be open in the same Matlab instance
    *   MMM is reinitialized, i.e., all models and existing visualization are deleted
	
``figures``
---------------------------------

Sets a general output format for saving graphics. 

.. code-block:: matlab

    figures format

Arguments
    *   ``format`` - figure output format, such as `png`, `pdf`, `bmp`, `jpeg`, `tiff`
Remarks
    *   default is `png` 
    *   figure format can also be specified in individual ``graphics`` commands
	
``getAlphaFold``
---------------------------------

Import an AlphaFold prediction via its UniProt identifier. 

.. code-block:: matlab

    getAlphaFold UniProtID

Arguments
    *   ``UniProtID`` - UniProt identifier, e.g. ``P61626``
Remarks
    *   note that not all proteins in UniProt have an AlphaFold prediction in the database
    *   there should be only one ``getAlphaFold`` line, or, alternatively, a ``getens``, an ``addpdb``, or an ``import`` line 
	
``getens``
---------------------------------

Input of an ensemble in MMMx ensemble list format. 

.. code-block:: matlab

    getens file

Arguments
    *   ``file`` - file name, extension .ens is appended if there is none
Remarks
    *   there should be only one ``getens`` line, or, alternatively, an ``addpdb`` line, an ``import``, or a ``getAlphaFold`` line 
    *   ``getens`` can also import from a '.zip' archive, as generated with the ``archive`` keyword of the EnsembleAnalysis keyword 
	
``getPED``
---------------------------------

Import an ensemble from the Protein Ensemble Database (PED). 

.. code-block:: matlab

    getPED PEDid

Arguments
    *   ``PEDid`` - PED identifier, e.g. ``PED00160.e001``, the part after the full stop is only needed if the entry contains several ensembles
Remarks
    *   one visualize block can only process only one ensemble, if the PED entry has several ensembles, ``.e00n`` is mandatory 
    *   do not combine with other keys that import an ensemble 
		
``getZenodo``
---------------------------------

Import an ensemble from Zenodo. 

.. code-block:: matlab

    getZenodo ZenodoID

Arguments
    *   ``ZenodoID`` - Zenodo identifier, e.g. ``6384003.hnRNPA1_unrestrained_raw_ensemble.zip``, the part after the full stop is the file name within the Zenodo entry
Remarks
    *   the Zenodo file can be a single PDB file or a ZIP archive or a (gzipped) TAR archive. 
    *   for a single PDB file, models are taken as conformers.
    *   if an archive contains several PDB files, they are taken as conformers
    *   an archive may additionally contain an MMMx ensemble (.ens) file. In this case, the ensemble is constructed according to this file
    *   a gzipped input file must contain only a single tar archive
		
``graphics``
---------------------------------

Request to save a graphic to a file or to copy it to the clipboard  

.. code-block:: matlab

    graphics [file [mode [view]]]
	
Arguments
    *   ``file`` - file name for the graphics file, must include extension if you need one 
    *   ``mode`` - graphics mode, such as `png`, `pdf`, `bmp`, `jpeg`, `tiff` 
    *   ``view`` - specification of viewing direction 

Remarks
    *   if there are no arguments, the current graphics is copied as a bitmap to the clipboard (Windows only)
    *   if the graphics mode is missing, it is specified by the ``figures`` keyword; if this is also missing, it is `png`
    *   ``view`` can be a Cartesian direction (`x`, `-x`, `y`, `-y`, `z`, `-z`) 
    *   alternatively, ``view`` can specifiy a viewing vector by three numbers, for instance `0.707 0.707 0` for halfway between `x` and `y`
    *   ``view`` can also be specified by six numbers; then, the final three numbers define the camera up direction, which must not coincide with the view direction
    *   use ``symmetry`` or ``bilayer`` in the ``prepare`` module for convenient coordinate transformations
    *   if you have a template with your preferred viewing orientation, use ``superimpose`` in the ``EnsembleAnalysis`` module for transformation
	
``import``
---------------------------------

Import an ensemble from PDB via its PDB identifier. This can also be used for loading a single PDB file. 

.. code-block:: matlab

    import pdbid

Arguments
    *   ``pdbid`` - PDB identifier, e.g. ``2LZM``, can also be a PDB file name, but then must have extension `.pdb`
Remarks
    *   if the PDB file has several conformers, all populations are set equal
    *   there should be only one ``import`` line, or, alternatively, a ``getens``, an ``addpdb``, or a ``getAlphaFold`` line
    *   ``getpdb`` is synonymous to ``import``	
	
``isosurface``
---------------------------------

Stand-alone isosurface visualization for density and property files (does not require MMM). This is a block key with options. 

.. code-block:: matlab

    isosurface density-file [property-file]
       option1 argument1 [argument2 argument3]
       ...
    .isosurface

Arguments
    *   ``density-file`` - name of the density file, use the EnsembleAnalysis module to generate one
    *   ``property-file`` - optional name of the property file for isosurface coloring, use the EnsembleAnalysis module to generate one
Available subkeys (options)
    *   ``colorscheme`` - property related color scheme, can be ``electrostatic`` (default), ``cation-pi``, or ``hydrophobic``
    *   ``level`` - fraction of total density included by the isosurface, defaults to 0.999, which is appropriate for ensemble pseudo-electron density
    *   ``camvec`` - vector pointing from isosurface to camera, three values, defaults to '1 0 0'
    *   ``camupvec`` - vector indicating the top of the camera, defaults to '0 1 0'
    *   ``limits`` - property level corresponding to extrema of the color scale, default depends on selected color scheme
    *   ``figname`` - figure name for saving, extension determines graphics format, defaults to 'isosurface.png', default extension is '.png'
    *   ``opaqueness`` - opaqueness of the isosurface, 0 is invisible, 1 is fully opaque, defaults to 1
Remarks
    *   if no property file is specified, the density isosurface is colored uniformly with SVG color 'gainsboro' 
    *   if a property file is specified, but no color scheme is specified, the Matlab default scheme 'parula' applies
    *   if a property file is specified, but no color scheme and no limits are specified, the limits of the color scale are the minimum and maximum property value 
	
``label``
---------------------------------

Generates and attaches spin label rotamers for later visualization.

.. code-block:: matlab

    label address type

Arguments
    *   ``address`` - MMM address of a residue (labelling site), e.g. `(A)131`
    *   ``type`` - label type, defaults to `mtsl`
Remarks
    *   the label is only generated, not shown, use ``show address label`` for visualizing the rotamer cloud

``normalize``
---------------------------------

Normalization mode for translation of populations to transparency. If `on` (default), 
the conformer with maximum population is completely opaque and opaqueness of other conformers is proportional
to the ratio of their population to the maximum population. If `off`, opaqueness equals population. 
The same normalization applies to coil radius in mode ``snake`` of keyword ``show``.

.. code-block:: matlab

    normalize mode

Arguments
    *   ``mode`` - can be `on` or `off`
    *   default is `on`

``script``
---------------------------------

Sets name of the MMM script file. Extension `.mmm` is added if there is none. 

.. code-block:: matlab

    script file

Arguments
    *   ``file`` - file name for the script file
Remarks
    *   default is MMMx.mmm
	
``show``
---------------------------------

Define graphics elements. All conformers in the ensemble have the same graphics elements; only transparency or width differs.

.. code-block:: matlab

    show address mode

Arguments
    *   ``address`` - MMM address of chains, residues, or atoms, use `(:)` for all chains
    *   ``mode`` - display mode, e.g., `ribbon`, such as ``graphics (:) ribbon``
Remarks
    *   all functionality of the ``show`` command of MMM is available
    *   in general, opacity (1-transparency) is proportional to population of conformers 
    *   an additional mode ``snake`` displays coils with radius proportional to population 
    *   in ``snake`` mode, all conformers are fully opaque, transparency is not used, this is faster 

