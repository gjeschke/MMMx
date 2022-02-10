.. MMMx documentation master file, created by
   sphinx-quickstart on Thu Jul 23 11:16:17 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

MMMx 1.0 Documentation
========================

*    MMMx is the successor program of MMM (Multiscale Modeling of Macromolecules)

*    it implements pipelines for integrative ensemble modeling of structures of proteins and protein complexes 

*    modeling is based on distance distribution information for determining the width of ensembles and for enhanced sampling in raw ensemble generation 

*    MMMx is an open-source Matlab package and is also available as compiled stand-alone Windows application

*    it consists of modules that can be combined to modeling pipelines in control files (MMMx scripts) 


.. image:: ../images/concept.png
   :width: 400
   :alt: module concept image
   
MMMx is available at `epr.ethz.ch <https://epr.ethz.ch/software.html>`_ and `at GitHub <https://github.com/gjeschke/MMMx>`_

.. toctree::
    :hidden:
    :caption: Installation
    :maxdepth: 2

    ./source_installation
    ./executable_installation


.. toctree::
    :hidden:
    :caption: Modules
    :maxdepth: 2

    ./experiment_design
    ./prepare
    ./rigi
    ./flex
    ./flex_RNA
    ./enm
    ./ensemble_fit
    ./ensemble_analysis
    ./locate
    ./yasara_refine
	
.. toctree::
    :hidden:
    :caption: Input/output concept
    :maxdepth: 2

    ./control_files
    ./log_files
    ./ensemble_list
    ./MMMx_addresses
	
.. toctree::
    :hidden:
    :caption: Spectroscopic labels
    :maxdepth: 2

    ./rotamer_concept
    ./distance_distributions

.. toctree::
    :hidden:
    :caption: Third-party software
    :maxdepth: 2

    ./third_party
    ./concept
	
.. toctree::
    :hidden:
    :caption: Examples
    :maxdepth: 2

    ./demo_ExperimentDesign
	

.. Indices and tables
.. ==================

.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`
