.. _control_files:

Control files
====================

Concept
---------------------------------

Modelling with MMMx can be performed either from the command line or from inside the GUI using the ``MMModel`` function.
``MMModel`` accepts a single parameter, which is the file name of a control file:

.. code-block:: matlab

    exceptions = MMModel(control_file)
    [exceptions,entity] = MMModel(control_file)

Parameters
    *   ``control_file`` - MMMx control file as described below

Returns
    *   ``exceptions`` - error message (cell array with MException objects for warnings or errors)
    *   ``entity``     - MMMx entity containing the model, may be populated by only some of the :ref:`modelling modules<modelling_modules>`

When MMModel is activated from the GUI, the control file can be selected with a file browser and exceptions are reported through the status text field of the GUI.

Any call to MMModel generates a :ref:`log file<log_files>`. The default file name is ``MMModel_yyyy_mm_dd_HH_MM.log``, where ``yyyy``, ``mm``, ``dd``, ``HH``, ``MM``
specify the year, month, day, hour, and minute of the starting time of ``MMModel``. A different log file name can be specified in the control file.

Further output files may be generated depending on the :ref:`modelling modules<modelling_modules>` that are called in the control file.

The control file specifies a modelling task by a sequence of module calls, each having its own run options and list of restraints. 
Where possible, MMMx checks whether restraints in a sequence of module calls are consistent.

Syntax
--------------

The syntax of MMMx control files is akin of the syntax of MMM restraint files, but MMM restraint files cannot be used without editing as MMMx control files.

The MMMx control file consists of blocks corresponding to individual :ref:`modelling modules<modelling_modules>`. Outside such blocks, keywords for input and output are supported.
Keywords outside of blocks have the form ``#`` *keyword*.

Comments are initiated by ``%`` and can appear anywhere. The rest of the line after a ``%`` character is ignored.

Module blocks are opened by ``!`` *module*, where *module* is the name of the :ref:`modelling module<modelling_modules>`. They are closed by ``.`` *module*. 
If a new block is opened or the file ends before the closing line of a block, ``MMModel`` proceeds, but raises a warning.

The syntax

``!`` *module* [*entity_1* [*entity_2*] ...] passes entities to the module, where *entity_1*, *entity_2* etc. are entity numbers. By default, a module requiring one entity will be passed the last entity that was used or generated.

The syntax  

``.`` *module* [*entity_1* [*entity_2*] ...] passes entities from the module to ``MMModel``, where *entity_1*, *entity_2* etc. are internal entity numbers.

The control file ends with a ``# end`` statement. If it is missing, MMMx raises a warning. All lines after an ``# end`` statement are ignored.

Each module can define its own keywords for restraint specification and run options. In principle, the same keyword can be used with different meanings in different modules, except for restraint keywords. 
However, such practice is discouraged.
Keywords inside blocks must be indented by at least two spaces. If a keyword initiates a block of restraint lines, the restraint lines must be indented by at least two spaces with respect to this keyword.
It is allowed and good practice to end such a block with ``.`` *keyword*.

Input and output specification
------------------------------

``# logfile`` *logname*  - name of the logfile, if present, first line of the control file

``# getpdb``  *pdbfile* [*number*] - fetches an entity from the PDB server or from a local file. *number* assigns an internal entity number (defaults to next available number)

``# putpdb``  *pdbfile* [*number*] - saves an entity to a local file. *number*  specifies the internal entity number (defaults to the last entity used or generated)

Restraint specification
-----------------------

Restraint format is the same for all modules.


