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

The control file can have any extension. If it has none and such a file without extension does not exist or is in the wrong format, ``.mxc`` is tried.

Any call to MMModel generates a :ref:`log file<log_files>`. The default file name is ``MMModel_yyyy_mm_dd_HH_MM.log``, where ``yyyy``, ``mm``, ``dd``, ``HH``, ``MM``
specify the year, month, day, hour, and minute of the starting time of ``MMModel``. A different log file name can be specified in the control file.

Further output files may be generated depending on the :ref:`modelling modules<modelling_modules>` that are called in the control file.

The control file specifies a modelling task by a sequence of module calls, each having its own run options and list of restraints. 
Where possible, MMMx checks whether restraints in a sequence of module calls are consistent.

Syntax
--------------

The syntax of MMMx control files is akin of the syntax of MMM restraint files, but MMM restraint files cannot be used without editing as MMMx control files.
Keywords are case-insensitive.

The MMMx control file consists of blocks corresponding to individual :ref:`modelling modules<modelling_modules>`. Outside such blocks, keywords for input and output are supported.
Keywords outside of blocks have the form ``#`` *keyword*.

Comments are initiated by ``%`` and can appear anywhere. The rest of the line after a ``%`` character is ignored.

Module blocks are opened by ``!`` *module*, where *module* is the name of the :ref:`modelling module<modelling_modules>`. They are closed by ``.`` *module*. 
If a new block is opened or the file ends before the closing line of a block, ``MMModel`` proceeds, but raises a warning.

The syntax

.. note::

  ``!``\ *module* [*entity_1* [*entity_2*] ...] 
  
passes entities to the module, where *entity_1*, *entity_2* etc. are entity numbers. By default, a module requiring one entity will be passed the last entity that was used or generated.

The syntax  

.. note::

  ``.``\ *module* [*entity_1* [*entity_2*] ...] 

passes entities from the module to ``MMModel``, where *entity_1*, *entity_2* etc. are internal entity numbers.

The control file ends with a ``# end`` statement. If it is missing, MMMx raises a warning. All lines after an ``# end`` statement are ignored.

Each module can define its own keywords for restraint specification and run options. In principle, the same keyword can be used with different meanings in different modules, except for restraint keywords. 
However, such practice is discouraged.
Keywords inside blocks must be indented by at least two spaces. If a keyword initiates a block of restraint lines, the restraint lines must be indented by at least two spaces with respect to this keyword.
It is allowed and good practice to end such a block with ``.`` *keyword*.

Input and output specification
------------------------------

.. note::

  ``#logfile`` *logname*  % name of the logfile

  ``#getpdb``  *pdbfile* [*number*] % fetch entity *number* from the PDB server or from a local file

  ``#putpdb``  *pdbfile* [*number*] % saves entity *number* to local file
  


Restraint specification
-----------------------

The following restraint formats apply to all modules. There can be multiple blocks of restraints of the same type within the same module section.
Please note that there may be additional module-specific restraint types that are explained in keyword specifications of the :ref:`modelling modules<modelling_modules>`.

**Distance distribution restraints**

Keyword: ``ddr``  Legacy keyword: ``deer``

The distance unit is **Angstroem**.

Syntax:

.. note:: 

  ``ddr`` *label_1* [*label_2*] [-permute]

      *site_1* *site_2* *<r1>* *fwhm1* [*pop1* *<r2>* *fwhm2* ...] ``% example for (multi)Gaussian restraints``
	
      *site_1* *site_2* ``-lb`` *lower_bound* ``-ub`` *upper_bound* ``% example for lower-bound/upper-bound restraints``
	
      *site_1* *site_2* ``@`` *distribution_data* ``% example for parameter-free distributions``
	
      *site_1* *site_2* ``-rice`` *<r1>* *fwhm1* [*pop1* *<r2>* *fwhm2* ...] ``% example for (multi)Rice restraints``

      *site_1* *site_2* ``-skew`` *<r>* *fwhm* *skewness* ``% example for skew Gauss restraints``
   
  ``.ddr``
 	
Each combination of label types at the two sites requires its own ``ddr`` block. If both sites are labeled with the same label, it is sufficient to specify it once.
Possible label types *label_1* and *label_2* correspond to existing :ref:`rotamer libraries<label_set>`. 
The syntax ``atom.<atname>`` specifies an atom, for instance, ``atom.CA`` the CA atom of the addressed site.

If option ``-permute`` is present, labels 1 and 2 are attached to sites 1 and 2 in both bossible permutations. 

The two labelled sites are specified by :ref:`MMMx residue addresses<MMMx_addresses>` *site_1* and *site_2* . 
These addresses must refer to either a residue in an entity passed to the module or to a residue generated by the module.

For a single Gaussian restraint, the mean value *<r1>* and full width at half maximum *fwhm1* need to be specified. 
For multi-Gaussian restraints, populations *pop1*, *pop2*, ... need to be specified, except for the last component (the populations add to 1).

Upper/lower bound restraints require the option specifiers ``-lb`` and ``ub``. These are hard restraints. Conformers are rejected, if the simulated distance is outside bounds.
Use distributions if you expect that some restraints may be violated.

Parameter-free distributions require an ASCII data file with name *distribution_data*.dat (extension may be included) that has at least two columns. The first column is the distance axis (Angstroem units) and the second column is the probability per distance bin.
We advise to provide lower and upper confidence limits for the bin probabilities in columns 3 and 4. White space between ``@`` and *distribution_data* is allowed, but not required.  
Bin probabilities are automatically normalized to unity sum.

(Multi)-Rice distributions require the option specifier ``-rice``. Otherwise they work line (multi)-Gaussian distributions.

Skew Gauss distribution restraints require the option specifier ``-skew``. Only a single skew distribution is supported.

Note that you can provide any parametrized distance distribution by first converting it to a binned distance distribution and saving it as *distribution_data*.dat.
However, for parametrized models that are supported, it is more transparent to specify the restraints by the parameters.

**Symmetry-related distance distribution restraints**

Symmetry constraints are specified at restraint level, because they may refer to only part of an entity.
In combination with distribution restraints, this choice allows for symmetry disorder.

Keyword: ``ddr_sym``  Legacy keyword: ``oligomer``

The distance unit is **Angstroem**.

Syntax:

.. note:: 

  ``ddr_sym`` *label* *n* [*ox* *oy* *oz* [*dx* *dy* *dz*]] [``-all``]

      *site* *<r1>* *fwhm1* [*pop1* *<r2>* *fwhm2* ...] ``% example for (multi)Gaussian restraints``
	
      *site* ``-lb`` *lower_bound* ``-ub`` *upper_bound* ``% example for lower-bound/upper-bound restraints``
	
      *site* ``@`` *distribution_data* ``% example for parameter-free distributions``
	
      *site* ``-rice`` *<r1>* *fwhm1* [*pop1* *<r2>* *fwhm2* ...] ``% example for (multi)Rice restraints``
 
      *site* ``-skew`` *<r>* *fwhm* *skewness* ``% example for skew Gauss restraints``
   
  ``.ddr_sym``
 	
Only one labeling site is specified. The second site is generated by a rotation by an angle 360°/*n* about a C\ :sub:`n` (*n*-fold) rotation axis.
In all-distances mode, specified by the option ``-all``, all *n(n-1)/2* pairs between the *n* label positions are considered. 
Otherwise, only the modal distance (side length of the *n* -gon) is considered (default). For *n* = 2 and *n* = 3, option ``-all`` has no effect.  

By default, the rotation axis is a *z* axis passing through the origin of the coordinate frame of the entity. 
It is possible to specify a point on the rotation axis by coordinates *ox*, *oy*, and *oz* 
as well as a direction of the rotation axis by coordinates *dx*, *dy*, and *dz*.

Apart from specification of only one site, restraint lines have the same syntax as in ``ddr`` blocks.


**Paramagnetic relaxation enhancement (PRE) restraints**

Keyword: ``pre``

The distance unit is **Angstroem**.

Syntax:

.. note:: 

  ``pre`` *label* *atom* [*taui* [*taur* [*taus*]]]

      *site_1* *site_2* *ratio* ``% example for ratio Ipara/Idia``
	
      *site_1* *site_2* ``-Gamma2`` *Gamma2* ``% example transverse relaxation enhancement``
	
  ``.pre``
 	
Possible spin label types *label* correspond to existing :ref:`rotamer libraries<label_set>`. If the entity has explicit protons, they can be specified by *atom*.
Otherwise, the heavy atom, for instance `N` for a backbone NH atom should be specified. MMMx will attempt to generate the proton position.

Correlation times can be provided as additional arguments. The default for the correlation time of internal motion (*taui*) is 250 ps, 
the one for global tumbling of the protein (*taur*) 3 ns, and the one for spin label relaxation (*taus*) 1 :math:`\mu`\ s.

The experimental restraints can be specified either by the intensity *ratio* between the paramagnetically and diamagnetically labelled sample or by the transverse relaxation enhancement rate *Gamma2*.


**SAXS restraints**

Small-angle x-ray scattering fits can be specified for the whole entity or for a subset of chains. 

Syntax:

.. note:: 

  ``saxs`` *saxs_data* [*sm*] [``-v3``] 

      *chain_1* [*chain_2* ...]
	
  ``.saxs``

Usually, SAXS restraints are specified by a single line, giving only the file name of the SAXS data and optionally a maximum scattering vector *sm* for fitting till *sm*.
SAXS fitting in MMMx uses `crysol`_ of the ATSAS package. Option ``-v3`` specifies that `crysol3`_ is used instead.
It is possible to specify only a subset of the chains of an enity for SAXS fitting. 
For this, ``saxs`` is used as a block restraint with a single additional line that specifies the included chains by :ref:`MMMx chain addresses<MMMx_addresses>`.

.. _crysol: https://www.embl-hamburg.de/biosaxs/manuals/crysol.html

.. _crysol3: https://www.embl-hamburg.de/biosaxs/manuals/crysol3.html

**SANS restraints**

Small-angle neutron scattering fits can be specified for the whole entity or for a subset of chains. 
If an experimental selection of part of the entity was made by contrast matching, 
it is better to specify the deuterium content in the buffer than to specify the selected chains. 

Syntax:

.. note:: 

  ``sans`` *sans_data* [*illres* [``*D2O*``]]

      *chain_1* [*chain_2* ...]
	
  ``.sans``

Usually, SANS restraints are specified by a single line, giving only the file name of the SANS data and optionally the name of a resolution file, 
a fraction of D\ :sub:`2`\ O in the solution. 
If the second argument is a number instead of a string, it is interpreted as D\ :sub:`2`\ O content and the resolution file is considered to be missing (a warning is raised).

SAXS fitting in MMMx uses `cryson`_ of the ATSAS package. 
It is possible to specify only a subset of the chains of an enity for SANS fitting. 
For this, ``sans`` is used as a block restraint with a single additional line that specifies the included chains by :ref:`MMMx chain addresses<MMMx_addresses>`.

.. _cryson: https://www.embl-hamburg.de/biosaxs/manuals/cryson.html

**Crosslink restraints**

Crosslink restraints can be specified as a fraction of potentially crosslinkable residue pairs that are sufficiently close to be actually crosslinked.
The maximum distance (in Angstroem) and fraction (0... 1) apply to all crosslinks in one block. 

Syntax:

.. note:: 

  ``crosslink`` *maxdist* *fraction* [*atom_a* [*atom_b*]]

      *site_1_a* *site_1_b*

      *site_2_a* *site_2_b*
	
      ...
	
  ``.crosslink``

The distance is measured between CA atoms, unless the atom types in sites a and b of the crosslink are specified. 
The ``crosslink`` line is followed by *x* lines specifying individual pairs of residues for which crosslinks were found.
The two linked sites are specified by :ref:`MMMx residue addresses<MMMx_addresses>`

A conformer is rejected if more than *fraction*\ ·\ *x* of the addressed atom pairs have a larger distance than *maxdist*.   

**Immersion depth restraints**

The depth of immersion of sites into a lipid bilayer can be specified with this restraint. 

Syntax:

.. note:: 

  ``depth`` *label* [*ox* *oy* *oz* [*dx* *dy* *dz*]]

      *site* *<r>* *fwhm* ``% example for Gaussian restraints``
	
      *site* ``-lb`` *lower_bound* ``-ub`` *upper_bound* ``% example for lower-bound/upper-bound restraints``
	
  ``.depth``

The depth is measured as a distance from the center plane of the bilayer, i.e., large values correspond to low immersion depth or even positions outside the bilayer.

Possible label types *label* correspond to existing :ref:`rotamer libraries<label_set>`. 
The syntax ``atom.<atname>`` specifies an atom, for instance, ``atom.CA`` the CA atom of the addressed site.

For a Gaussian restraint, the mean value *<r>* and full width at half maximum *fwhm* need to be specified.  
Upper/lower bound restraints require the option specifiers ``-lb`` and ``ub``. These are hard restraints. 
Conformers are rejected, if the simulated distance is outside bounds. Use distributions if you expect that some restraints may be violated.

By default, the bilayer normal is assumed to be the *z* axis and the center plane is assumed to pass through *z* = 0 of the coordinate frame of the entity. 
It is possible to specify a point on the center plane by coordinates *ox*, *oy*, and *oz* and the direction of the bilayer normal by coordinates *dx*, *dy*, and *dz*.

**Secondary structure and cis peptide propensities**

The keywords ``alpha``, ``beta``, ``polypro``, and ``cis`` allow to specify propensities at a residue to adopt :math:`{\alpha}`\ -helix, :math:`{\beta}`\ -strand, polyproline-helix, or cis-peptide backbone torsion angles. 
The following example is for :math:`{\alpha}`\ -helix propensities.

Syntax:

.. note:: 

  ``alpha``

      *site* *propensity* ``% example for a single site``
	
	  *site_1*\ ``-``\ *site_2* *propensity* ``% example for a range of residues``
	
  ``.alpha``

The sites *site*, &site_1*, and *site_2* are specified by :ref:`MMMx residue addresses<MMMx_addresses>` and *propensity*  is a value between 0 and 1.
Use the range syntax with propensity 1 to strictly enforce secondary structure for a certain section of residues. 

Specifying propensities instead of physical ensemble mean restraints related to them (e.g. NMR chemical shifts and residual dipolar couplings) 
is prefereable in ensemble building, as it allows to adapt backbone torsions statistics, which in turn improves sampling of suitable conformations.

In ensemble fitting, it is advisable to specify restraints as close as possible to primary experimental data.

**Rigid bodies**

In general, a rigid body can comprise one or more sections of one or more macromolecular chains, as well as cofactors or other ligands.
Builder modules, such as ``Rigi``, may require that the complete chain behaves as a rigid body. Template entities have to be prepared to ensure this.

Syntax:

.. note:: 

  ``rigid`` *section_1* [*section_2* [...]]

      *refsite_1* *reflabel_1* ``% origin``
	
      *refsite_2* *reflabel_2* ``% point on *x* axis``
	
      *refsite_3* *reflabel_3* ``% point in *xy* plane``
	
  ``.rigid``

Rigid bodies are internally numbered in the sequence of the corresponding ``rigid`` blocks. This numbering is internal to a module.

The section specifiers *section_1*, *section_2*, ... are :ref:`MMMx chain or residue-range addresses<MMMx_addresses>`, such as ``(B)`` or ``(C)58-123``.

The three reference sites *refsite_1*, *refsite_2*. and *refsite_3* are obligatory and should not be situated on a line. 
They specify a local frame and can be used for computing rigid-body arrangements by distance geometry.  

The labels *reflabel_1*, *reflabel_2*, and *reflabel_3* either correspond to existing :ref:`rotamer libraries<label_set>` or 
have the syntax ``atom.<atname>``. The latter syntax specifies an atom, for instance, ``atom.CA`` the CA atom of the reference site.