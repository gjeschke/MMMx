function [taur,info] = rotational_correlation_time(entity,address,options)
%
% ROTATIONAL_CORRELATION_TIME   Returns rotational correlation time of a
%                               chain or complex
%
%   [taur,info] = ROTATIONAL_CORRELATION_TIME(entity)
%   Given an entity, a vector of rotational correlation times is returned
%   for all conformers
%
%   [taur,info] = ROTATIONAL_CORRELATION_TIME(entity,address)
%   Returns the rotational correlation time for an addressed object in an
%   entity
%
%   [taur,info] = ROTATIONAL_CORRELATION_TIME(entity,address,options)
%   Allows to specify options, with empty address, a vector of rotational
%   correlation times is returned for all conformers in the entity
%
% INPUT
% entity       entity in an MMMx format, must be provided
% address      MMMx address, 'selected' refers to the current selection
% options      struct, the following options can be specified
%              .viscosity   viscosity in cPoise, defaults to 0.8882 cPoise,
%                           (pH 7.4 phosphate buffer at 25C)
%              .T           temperature, defaults to 25C, can be given in
%                           ?C or K, values above 100 default to K units
%              .specV       partial specific volume, defaults to 0.702
%                           g/cm^3, used only with HYDROPRO
%              .density     solvent density, defaults to 1 g/cm^3, used
%                           only with HYDROPRO
%              .HYDROPRO    Boolean flag, if true, HYDROPRO is used (if on
%                           Matlab path), defaults to false
%              .N           number of residues, must be provided if
%                           HYDROPRO is not true
%
% OUTPUT
% taur         rotational correlation time or vector thereof
% info         struct, type of computation and further results
%              .HYDROPRO    Boolean flag, if true, HYDROPRO was used
%              .Rh          hydrodynamic radius 
%              .Tr          harmonic mean relaxation time
%
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2020: Gunnar Jeschke

%set empty output
taur = [];
info.HYDROPRO = false;
info.Rh = [];
info.error = '';

if ~exist('options','var') || ~isfield(options,'HYDROPRO') || isempty(options.HYDROPRO)
    options.HYDROPRO = false;
end

if options.HYDROPRO
    options.executable = which_third_party_module('hydropro10');
    if ~isempty(options.executable)
        info.HYDROPRO = true;
    end
end

if ~info.HYDROPRO 
    if ~isfield(options,'N') || isempty(options.N)
        info.error = 'Number of residues is missing';
        return
    end
end

if ~isfield(options,'viscosity') || isempty(options.viscosity)
    options.viscosity = 0.8882;
end

if ~isfield(options,'T') || isempty(options.T)
    options.T = 273.15+25;
end

% correct Celsius temperatures to Kelvin
if options.T <= 100
    options.T = options.T + 273.15;
end

if ~isfield(options,'density') || isempty(options.density)
    options.density = 1.0;
end

if ~isfield(options,'specV') || isempty(options.specV)
    options.specV = 0.702;
end



overwrite = true;
if ~exist('address','var') || isempty(address) % all conformers of the entity
    C = length(entity.populations);
    taur = zeros(1,C);
    info.Rh = zeros(1,C);
    for c = 1:C
        entity = select(entity,sprintf('{%i}(*).*',c),overwrite);
        [taurc,Rhc,Tr] = get_taurc(entity,options,info);
        taur(c) = taurc;
        info.Rh(c) = Rhc;
        info.Tr(c) = Tr;
    end
else
    entity = select(entity,address,overwrite);
    [taur,Rhc,Tr] = get_taurc(entity,options,info);
    info.Rh = Rhc;
    info.Tr = Tr;
end

% exceptions = put_pdb(entity,fname,options)

function [taur,Rhc,Tr] = get_taurc(entity,options,info)

% empty output
taur = [];
Rhc = [];

kB = 1.38064852e-23; % Boltzmann constant

pse= load('element_attributes.mat');

if info.HYDROPRO
    save_options.selected = true;
    [HP_path,~] = fileparts(options.executable);
    my_path = pwd;
    cd(HP_path);
    [~,~,elements] = get_coor(entity,'selected');
    MW = sum(pse.pse.mass(elements));
    put_pdb(entity,'mmmx.pdb',save_options);
    fid = fopen('hydropro.dat','wt');
    fprintf(fid,'MMMx-rotcorr\n'); 
    fprintf(fid,'MMMx_HYDROPRO\n'); 
    fprintf(fid,'%s\n','mmmx.pdb');
    fprintf(fid,'1\n');
    fprintf(fid,'2.9,\n');
    fprintf(fid,'6,\n');
    fprintf(fid,'1.0,\n');
    fprintf(fid,'2.0,\n');
    fortran = sprintf('%6.1f,',round(options.T-273.15));
    fprintf(fid,'%s\n',strtrim(fortran));
    fortran = sprintf('%8.4f,\n',options.viscosity/100);
    fprintf(fid,'%s\n',strtrim(fortran));
    fortran = sprintf('%8.1f,',MW);
    fprintf(fid,'%s\n',strtrim(fortran));
    fortran = sprintf('%6.3f,',options.specV);
    fprintf(fid,'%s\n',strtrim(fortran));
    fortran = sprintf('%6.3f,',options.density);
    fprintf(fid,'%s\n',strtrim(fortran));
    fprintf(fid,'-1\n');
    fprintf(fid,'-1\n');
    fprintf(fid,'0,\n');
    fprintf(fid,'1\n');
    fprintf(fid,'*\n');
    fclose(fid);
    s = dos(options.executable);
    if s == 0
        fid = fopen('MMMx_HYDROPRO-res.txt','rt');
        while 1
            tline = fgetl(fid);
            if ~ischar(tline)
                break
            end
            args = split(tline,':');
            disp(args{1});
            if strcmp(strtrim(args{1}),'Harmonic mean (correlation) time')
                nargs = split(args{2});
                Tr = str2double(nargs{2});
            end
            if strcmp(strtrim(args{1}),'Rotation (h)')
                nargs = split(args{2});
                Rhc = 1e8*str2double(nargs{2});
            end
        end
        eta = 0.001*options.viscosity; % cPoise to Pa s
        Rh = Rhc*1e-10; % Angstroem to m
        taur = 4*pi*eta*Rh^3/(3*kB*options.T);
        fclose(fid);
    end
    cd(my_path);
else
    coor = get_coor(entity,'selected');
    Rhc = hydrodynamic_radius(coor,options.N); % hydrodynamic radius
    % use Stokes-Einstein equation
    eta = 0.001*options.viscosity; % cPoise to Pa s
    Rh = Rhc*1e-10; % Angstroem to m
    taur = 4*pi*eta*Rh^3/(3*kB*options.T);
end