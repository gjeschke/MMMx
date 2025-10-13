function mk_OpenMM_minimize_script(fname,filelist,pH)
% mk_OpenMM_minimize_script(fname,filelist)
%
%    Generates a Python script for minimizing individual PDB structures in
%    an ensemble by OpenMM. The protein structures are solvated in an
%    aqueous solution with 150 mM NaCl with pH 7.0 (default) and
%    neutralized by counterions
%
% Input:
%
% fname     file name for the script, extension .py is appended if
%           extension is missing
% filelist  array of structures with single field (to be consistent
%           with dir)
%               .name file name
% pH        pH value of the solvent, defaults to 7.0
%
% G. Jeschke, 13.10.2025

if ~exist('pH','var') || isempty(pH)
    pH = 7.0;
end

[path,bname,ext] = fileparts(fname);
if isempty(ext)
    fname = fullfile(path,strcat(bname,'.py'));
end

fid = fopen(fname,'w');

% write preamble that loads the required OpenMM stuff 
fprintf(fid,'from openmm.app import *\n');
fprintf(fid,'from openmm import *\n');
fprintf(fid,'from openmm.unit import *\n\n');

% define forcefield
fprintf(fid,'forcefield = ForceField(''amber14-all.xml'', ''amber14/tip3pfb.xml'')\n');

% write the minimization parts for all input structures
for c = 1:length(filelist)
    fprintf(fid,'# Minimization of %s\n\n',filelist(c).name);
    [path,basname,~] = fileparts(filelist(c).name);
    outname = fullfile(path,sprintf('%s_optimized.pdb',basname)); 
    % Load this PDB file
    fprintf(fid,'pdb = PDBFile(''%s'')\n',filelist(c).name);
    % create modeller and add hydrogens
    fprintf(fid,'modeller = Modeller(pdb.topology, pdb.positions)\n');
    fprintf(fid,'modeller.addHydrogens(forcefield, pH=%4.1f)\n',pH);
    % solvate and add ions (neutralize and set [NaCl] = 0.15 mol/L)
    fprintf(fid,'modeller.addSolvent(forcefield,\n');
    fprintf(fid,'                model=''tip3p'',\n');
    fprintf(fid,'                padding=1.0*nanometer,\n');
    fprintf(fid,'                ionicStrength=0.15*molar,\n');
    fprintf(fid,'                neutralize=True,\n');
    fprintf(fid,'                positiveIon=''Na+'',\n');
    fprintf(fid,'                negativeIon=''Cl-'')\n');
    % System setup for simulation/minimization
    fprintf(fid,'system = forcefield.createSystem(modeller.topology,\n');
    fprintf(fid,'                             nonbondedMethod=PME,\n');
    fprintf(fid,'                             nonbondedCutoff=1.0*nanometer,\n');
    fprintf(fid,'                             constraints=HBonds)\n');
    fprintf(fid,'integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)\n');
    fprintf(fid,'simulation = Simulation(modeller.topology, system, integrator)\n');
    fprintf(fid,'simulation.context.setPositions(modeller.positions)\n');
    % energy minimization (optimization)
    fprintf(fid,'simulation.minimizeEnergy()\n');
    fprintf(fid,'state = simulation.context.getState(getPositions=True)\n');
    fprintf(fid,'positions = state.getPositions()\n');
    % create a Modeller from the full topology and minimized positions
    fprintf(fid,'modeller = Modeller(simulation.topology, positions)\n');
    % remove water molecules
    fprintf(fid,'modeller.deleteWater()\n');
    % remove ions
    fprintf(fid,'ions_to_remove = [res for res in modeller.topology.residues() if res.name in (''NA'', ''CL'')]\n');
    fprintf(fid,'modeller.delete(ions_to_remove)\n');
    % write the minimized protein-only structure to PDB
    fprintf(fid,'with open(''%s'', ''w'') as f:\n',outname);
    fprintf(fid,'     PDBFile.writeFile(modeller.topology, modeller.positions, f)\n\n');
end

fclose(fid);
