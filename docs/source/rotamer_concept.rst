.. _rotamer_concept:

Rotamer libraries
====================

Concept
---------------------------------

As previously established in MMM, site-directed labeling is modelled in MMMx based on pre-computed rotamer libraries of the label side groups.
Such libraries can be generated with classical atomic force fields via molecular dynamics simulations or Monte Carlo sampling of conformation space.
Briefly, the label side group is represented by a moderate number, typically between 100 and 10000, of rotameric states and their associated 
reference populations for a state that is non-interacting with the protein. Upon attachment to a labelling site, interaction with the protein is computed 
for a static (ensemble) structure of the latter assuming only non-bonded interaction via a Lennard-Jones potential parametrized in the Universal Force Field (UFF). 
Flexibility of the protein as well as of the label beyond rotameric states is roughly modelled by an *f*-factor (`0 \le f \le 1`) that scales van-der-Waals radii of atoms in the underlying force field. 
A few libraries use an enhanced attraction term, i.e., the attractive part of the Lennard-Jones potential is multiplied by the *f*-factor and an enhancement factor `e_\mathrm{attract}`. 
The default *f*-factor is 0.5 and the default enhancement factor is 1.0.

Attachment results in rotation and translation of the rotamer coordinates and in reweighting of populations, assuming a Boltzmann distribution at 298 K. 
Distributions of properties can then be computed by population-weighted averaging.

A number of features in MMM were used only in method development and have not been in use anymore for years. These features are deprecated in MMMx.
The concept of a "library temperature" has been deprecated as well, as libraries computed for the glass transition temperature of the medium have been found to perform worse than those computed for ambient temperature (298 K).

Implementation
--------------

The labeling concept is implemented by providing a set of rotamer libraries, a package for generating additional libraries for new labels, 
and by a function ``get_label`` for retrieving attributes of labels, such as distributions of label position and orientation.

Labelling is dynamic and virtual. Dynamic labeling means that a label *L* at residue *R* is computed at the time when the user of an entity first tries 
to retrieve attributes of this label at this particular site. The label attributes are then stored at residue level, whereas no atomic coordinates are 
generated at entity level (virtual labeling). Explicit coordinates of all atoms are generated only for clash tests, for saving label rotamers to a PDB file,
for visualizing the label in all-atom graphics, or when requested. Once computed, the explicit coordinates are stored in the atom coordinate array  of the
entity, while their indices are stored in the ``labels`` field of the labelled residue.

The following chemical modifications are implicit (d on need to be performed before labeling):

* mutation of any amino acid to cysteine or an unnatural amino acid for labeling

* conversion of uracil to 4-thiouracil

* conversion of phosphates to thiophosphates

.. _get_label:

The ``get_label`` function
----------------------------- 

Use this function if you want to label a residue with a particular label or want to retrieve attributes of a previously computed label.

.. code-block:: matlab

    argsout = get_label(entity,label,attributes)
    [argsout,exceptions] = get_label(entity,label,attributes)
    argsout = get_label(entity,label,attributes,address)
    [argsout,exceptions] = get_label(entity,label,attributes,address)


Parameters
    *   ``entity`` - entity in MMMx:atomic format
    *   ``label`` - label name
    *   ``attributes`` - either a string (see table below) or a cell array of strings
    *   ``address`` - MMMx residue address
Returns
    *   ``argsout`` - output arguments (*M*-element cell array for a single attribute)
    *   ``exceptions`` - error message (1-element cell array)
	
Either a rotamer library for ``label`` must be registered with MMMx (see below) or the label name must be ``atom.<atname>``, where ``<atname>`` is an existing atom name for this residue.

Requesting several attributes for the same site (residue) simultaneously can lead to a significant speedup in large entitities.
The gain in speed may be very large if the labels were already precomputed. To do this, arrange all attributes in a cell of strings and split the output cell array:

.. code-block:: matlab

    site = sprintf('{%i}%s',conformer_number,residue);
    [argsout,entity,exceptions] = get_label(entity,label,{'positions','populations'},site);
    positions = argsout{1}{1};
    populations = argsout{2}{1};
	
**Attributes**
	
====================== =============================================== =================
Variable               Explanation                                     Type   
====================== =============================================== =================
``info``               rotamer library information                     struct
``info.tlc``           MMMx internal three-letter code                 string
``info.class``         label class, for instance ``nitroxide``         string
``info.smiles``        SMILES string representing label structure      string
``info.f_factor``      f-factor and attraction enhancement factor      (1,2) double
``info.ref_pop``       reference populations for *R* rotamers          (*R*,1) double
``info.site``          site address                                    string
``part_fun``           attachment partitition function                 double
``potentials``         attachment potentials (J/mol)                   (*R*,1) double
``torsion``            *T* sidechain torsions for *R* rotamers         (*R,T*) double
``orientations``       Euler angles of molecular frame                 (*R*,3) double
``populations``        populations for *R* rotamers                    (*R*,1) double
``positions``          positions for *R* rotamers                      (*R*,3) double
``numbers``            numbers of the rotamers in the library          (*R*,1) int
``coor``               full atom coordinates of *R* rotamers           (*R*,1) cell
``affine``             affine matrix that transforms from the standard (4,4) double
                       frame to the site frame
====================== =============================================== =================
 
Note that SMILES strings for nitroxides tend to be interpreted as the corresponding hydroxylamines by some programs, notably by ChemDraw.
For libraries that contain several stereoisomers, the SMILES string refers to only one of them.

The attribute ``orientations`` can be used for simulating orientation selection in pulsed dipolar spectroscopy for spin labels
or `\kappa` averaging for FRET chromophores. In order to compute unit vectors along the `x,y,z` axes of the label molecular frame
in the entity frame (direction cosine matrix ``DCM``, the unit vectors are matrix rows) for rotamer `r=1`, use the code (example):

.. code-block:: matlab
     
   entity = get_pdb('2lzm'); % load structure of T4 Lysozyme with ID 2lzm from PDB server 
   [argout,exceptions] = get_label(entity,'mtsl','orientations','(A)131'); % get MTSL label orientation at residue 13
   orientations = argout{1}; % extract from cell output
   r = 1;
   DCM = Euler2DCM(orientations(r,:)); % compute direction cosine matrix

Set of libraries
-----------------

MMMx uses a single type of rotamer libraries, built by hierarchical clustering of Monte-Carlo generated conformer ensembles.
Ensemble generation assumes torsion potentials and the non-bonded interaction potential of UFF. 

For the following labels, rotamer libraries are provided with MMMx:

.. _label_set:

=======  ===============================  ==============  ================ ========= ====================
Label    Synonyms                         Class           Attachment       Rotamers  `e_\mathrm{attract}`
=======  ===============================  ==============  ================ ========= ====================
``R1A``  ``mtsl``, ``mtssl``              nitroxide       cysteine         216       1
``R7A``  ``br-mtsl``, ``br-mtssl``        nitroxide       cysteine         216       1
``V1A``  ``v1``                           nitroxide       cysteine         72        1
``IA1``  ``ia-proxyl``                    nitroxide       cysteine         108       1
``BAS``  ``basl``                         nitroxide       cysteine         240       1
``MA1``  ``ma-proxyl``                    nitroxide       cysteine         108       2.0
``DZD``  ``dzd``                          nitroxide       cysteine         216       1
``DZC``  ``dzc``                          nitroxide       cysteine         216       1
``GDI``  ``iag``                          nitroxide       cysteine         2461      2.0
``GDM``  ``mag``                          nitroxide       cysteine         1367      2.0
``TUP``  ``iap-4tu``, ``iap-thiouracil``  nitroxide       4-thiouracil     72        1                          
``TUM``  ``mts-4tu``, ``mts-thiouracil``  nitroxide       4-thiouracil     192       1
``RTT``  ``r5-tpt``                       nitroxide       5'-thiophosphate 576       1
``R5T``  ``r5-tp``                        nitroxide       thiophosphate    360       1
``R5P``  ``r5p``                          nitroxide       5'-thiophosphate 2048      1
``R3P``  ``r3p``                          nitroxide       3'-thiophosphate 512       1
``K1H``  ``HF-K1``                        nitroxide       unnatural aa     288       1
``HO4``  ``HO4451``                       nitroxide       unnatural aa     1024      1
``NC1``  ``cNox@Tyr``                     nitroxide       tyrosine         128       1
``NX1``  ``lNox@Tyr``                     nitroxide       tyrosine         256       1
``CNR``  ``CNC-NO``                       nitroxide       cofactor         144       1
``GMO``  ``dota-gd``                      gadolinium      cysteine         648       1
``GTO``  ``dtpa-gd``                      gadolinium      cysteine         2430      1
``GBH``  ``br-py-do3a``                   gadolinium      cysteine         180       1
``GBM``  ``br-py-do3ma``                  gadolinium      cysteine         180       1
``M8D``  ``m8-dota-gd``                   gadolinium      cysteine         1944      1
``GPM``  ``gpymi-MTA``                    gadolinium      cysteine         432       1
``TMT``  ``tormyshev-trityl``             trityl          cysteine         3888      1
``HCU``  ``dHis-Cu``                      histidine       any amino acid   12        1
=======  ===============================  ==============  ================ ========= ====================

Label names (three-letter codes) and synonyms are case-insensitive. 
Note that gadolinium labels are sufficiently good approximations for other lanthanide labels with the same ligand,
for instance, for pseudo-contact shift (PCS) and paramagnetic relaxation enhncement (PRE) computations. 

-----------------------------

Rotamer library format
----------------------

Rotamer libraries are stored in a binary Matlab format as a ``struct`` variable ``rot_lib``. The fields are defined as follows:

	
======================= =============================================== ================================
Field                   Content                                         Type   
======================= =============================================== ================================
``tlc``                 three-letter code                               string
``synonyms``            *S* synonyms for the label name                 (1,*S*) cell string
``SMILES``              SMILES string defining the structure            string
``rotamers``            information on *R* rotamers (index *r*)         (1,*R*) struct  
``rotamers(r).coor``    Cartesian coordinates of ``A`` atoms/rotamer    (*A*,3) double
``rotamers(k).torsion`` values of *T* torsion angles `\chi_t` for       (*R,T*) double
                        *R* rotamers
``elements``            atomic numbers for *A* atoms                    (*A*,1) uint8
``populations``         *R* populations for the non-attached rotamers   (*R*,1) double
``position``            *P* atom numbers and densities that define the  (*P*,2) double
                        label position
``attachment``          structure element to which the label can be     string
                        attached, for instance ``peptide``
``side_chain``          first side chain atom numer, for instance ``9`` ìnt
                        for a ``CB`` atom in position 9
``f_factor``            "forgive" factor and attraction enhancement     (1,2) double
``atom_tags``           *A* atom names                                  (1,*A*) cell string
``std_frame``           atoms that define the standard frame for
                        attachment: origin, atom on *x* axis, atom in   (1,3) int
                        *xy* plane
``std_frame_atoms``     atom types of the standard frame                (1,3) cell string 
``mol_frame``           atoms that define the label molecular frame:    (1,3) int
                        origin, atom on *x* axis, atom in *xy* plane
``mol_frame_atoms``     atom types of the standard frame                (1,3) cell string 
``class``               label class, for instance ``nitroxide``         string
``chi_def``             definition of *T* torsion angles `\chi_t`       (*T*,4) int
``connect``             bonding information for up to *B* bonds for *A* (*A,B*) int
                        atoms
``attach_forcefield``   force field for protein attachment, usually     string
                        ``UFF_Towhee``
``B_factors``           pseudo-temperature factors for *N* atoms in     (*N,R*) double
                        *R* rotamers
``method``              method for library generation. for instance     string
                        ``MMM_Monte_Carlo``
``gen_forcefield``      force field used in library generation          string
                        for instance ``UFF_Towhee``		
``types``               *A* atom type numbers for the used force field  (*A*,1) uint16 						
``solvation``           solvation assumed in library generation,        string
                        usually ``none``						
``prerun``              number of trials in a prerun of the MMMx native int
                        UFF Monte Carlo rotamer library generator
``suppress_H``          true if hydrogen atoms were neglected in        boolean
                        library generation, not recommended
``threshold``           thresholds for confermer acceptance in the      double
                        MMMx native generator
``min_strain``          minimum strain energy (kcal/mol) encountered    double
``maxdist``             maximum distance of the label position from     double
                        the backbone (CA atom or origin of attachment)
                        frame					   
``color``               RGB color triplet (fraction) for display        (1,3) double
``radius``              sphere radius (Angstroem) corresponding to      double
                        100% rotamer population for display
``ff``                  force field parameters for attachment           struct
``ff.LJ_r``             van-der-Waals radii for attachment indexed by   (1,103) double
                        atomic number
``ff.LJ_D``             Lennard-Jones potentials for attachment indexed (1,103) double
                        by atomic number
``ff.types``            atom type tags for force field                  string
======================= =============================================== ================================ 

The pseudo-temperature factors relate to the variation of atom positions within the cluster of conformers that was projected onto a single rotamer. 


