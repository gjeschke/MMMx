.. MMMx documentation master file, created by
   sphinx-quickstart on Thu Jul 23 11:16:17 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

MMMx 0.0.0 Documentation
========================

.. warning:: MMMx is at an early design and prototyping stage. Version 0.0.0 is not functional.

*    MMMx is the successor program of MMM (Multiscale Modeling of Macromolecules)

*    Graphical user interface is separated from computation modules

*    Methods can be implemented in any language. They communicate by a few defined data exchange formats

.. toctree::
    :hidden:
    :caption: Ensemble structure
    :maxdepth: 2

    ./functional_entity
    ./structure_representation
    ./MMMx_addresses
    ./MMMx_atomic
    ./MMMx_RigiFlex
	
.. toctree::
    :hidden:
    :caption: Spectroscopic labels
    :maxdepth: 2

    ./rotamer_concept
    ./distance_distributions

.. toctree::
    :hidden:
    :caption: Programmer's Guide
    :maxdepth: 2

    ./object_access
    ./object_modification
    ./property_functions
    ./ChimeraX
    ./python
    ./dependencies

.. toctree::
    :hidden:
    :caption: Architecture
    :maxdepth: 2

    ./concept
    ./data
    ./docs


.. Indices and tables
.. ==================

.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`
