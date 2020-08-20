.. _log_files:

.. highlight:: rst

Log files
====================

Concept
---------------------------------

Modelling in MMMx or use of the MMMx GUI automatically generate log files. These log files adopt the ReStructuredText markup text format, 
as it is easily readable by humans and at the same time can be easily parsed and converted into HTML website documents. 

Default log file names contain time information on their original creation, such as ``MMModel_yyyy_mm_dd_HH_MM.log`` or ``MMMxGUI_yyyy_mm_dd_HH_MM.log``, where ``yyyy``, ``mm``, ``dd``, ``HH``, ``MM``
specify the year, month, day, hour, and minute.

Log files keep the full information on a GUI session or modelling run. They can be processed into files at lower verbose level.

Conventions
--------------

Commands (::), module calls (!), and calls of external programs (>), such as ChimeraX, are echoed together with their arguments::

  :: `command` [`arg1`  [`arg2` ...]]

  ! `module` [`arg1` [`arg2` ...]] 

  > ChimeraX `command` [`arg1` [`arg2` ...]] 

Such lines open an indented block that contains all output generated during the task. Once the task is completed, the block is closed::

  .: `command`

  . `module` 

  < ChimeraX `command`

If a task generates a further call, it opens a new subblock. Example::

  ! EnsembleFit hnRNPA1_LCD.mcx

    Loaded model hnRNPA1_m1.pdb
  
    Processed 22 distance distribution restraints
  
    > crysol A1_minus_buffer_concnorm.dat
  
             chi^2 = 1.125
			   
    < crysol

    ...
  
  . EnsembleFit
