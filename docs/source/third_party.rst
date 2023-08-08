.. _third_party:

Third-party software
==========================

MMMx uses third-party software for tasks where well-tested solutions are already available in the form of executables. The executables must be on the Matlab path.  

ATSAS
----------------------------------

`ATSAS <https://www.embl-hamburg.de/biosaxs/software.html>`_ is used for simulating and fitting small-angle scattering (SAXS and SANS)
curves in module :ref:`EnsembleFit<ensemble_fit>`. The programs ``crysol`` and ``cryson`` from this package are accessed.

DSSP
----------------------------------

If `DSSP <https://swift.cmbi.umcn.nl/gv/dssp/HTML/distrib.html>`_ is present, MMMx can perform DSSP analysis upon loading 
PDB files. In version 1.0, this is of interest only for programmatic access. None of the current modules relies on DSSP secondary
structure information.

HYDROPRO
----------------------------------

`HYDROPRO <http://leonardo.inf.um.es/macromol/programs/hydropro/hydropro.htm>`_ can be used for computing
hydrodynamic radii as a pre-requisite for predicting paramagnetic relaxation enhancements (PREs) in module 
:ref:`EnsembleFit<ensemble_fit>`. This is rather slow and usually, a simpler estimate incorprated into MMMx suffices.

MSMS
----------------------------------

`MSMS <http://mgl.scripps.edu/people/sanner/html/msms_home.html>`_ is used for computation of solvent-accessible 
surfaces in lipid-bilayer optimization in  module :ref:`Prepare<prepare>`.

MUSCLE
----------------------------------

`MUSCLE <http://www.drive5.com/muscle/downloads.htm>`_ is used for sequence alignment in superposition of structures
in module :ref:`Prepare<prepare>`.

SCWRL4
---------------------------------

`SCWRL4 <http://dunbrack.fccc.edu/SCWRL3.php/>`_ is used for adding or modifying amino acid sidechains in module
:ref:`Flex<flex>` for loop modeling, in module :ref:`Prepare<prepare>` for mutations, sidechain repair, and sidechain repacking,
and in module :ref:`YasaraRefine<yasara_refine>` for optional repacking before calling YASARA in order to repair serious clashes
that may upset YASARA. Download from the Dunbrack lab requires a license application for non-profit users. On Linux systems, the path to SCWRL4 must be specified explicitly.
This is realized in the :ref:`Flex<flex>` module by the ``scwrl4`` keyword.

YASARA
----------------------------------

`YASARA structure <http://www.yasara.org/>`_ is the core of module :ref:`YasaraRefine<yasara_refine>`.
YASARA structure requires a paid license. In some cases, it is advisable to perform such refinement on
the raw ensemble before using module :ref:`EnsembleFit<ensemble_fit>`. As an alternative, the raw ensemble, which is a set 
of PDB files can be refined with external software. In many cases it is sufficient to refine the fitted ensemble.