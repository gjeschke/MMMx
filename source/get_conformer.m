function [argout,exceptions] = get_conformer(entity,attribute,conformer)
%
% GET_CONFORMER Retrieves attributes of a conformer
%
%   [argout,exceptions] = GET_CONFORMER(entity,attribute)
%   Provides attribute values and possibly exceptions for the first  
%   conformer of an entity
%
%   [argout,exceptions] = GET_CONFORMER(entity,attribute,conformer)
%   Provides attribute values and possibly exceptions for one conformer 
%   selected by address
%
% INPUT
% entity       entity in an MMMx format, must be provided
% address      MMMx address, 'selected' refers to the current selection
%              defaults to 'selected'
% attribute    chain attribute to be retrieved, string, defaults to
%              'info'
%              attribute   output                             Type
%              ------------------------------------------------------------
%              dssp        DSSP information for C chains      (1,C) struct
%              dssp(c)     .chain       chain identifier      string
%                          .sequence    aa sequence           string
%                          .secondary   secondary structure   string
%                          .sheets      sheet information     (R,2) double
%                          .bp          bridge partners       (R,2) double
%                          .NHO         N-H->O H bonds        (1,R) struct 
%                                       .hp   partner         int
%                                       .en   energy          double
%                          .OHN         O->N-H H bonds        (1,R) struct 
%                                       .hp   partner         int
%                                       .en   energy          double
%                          .tco         C=O cosine r,r+1      double
%                          .kappa       CA angle r-2,r,r+2    double
%                          .alpha       CA torsion r-1... r+2 double
%                          .phi         backbone dihedral phi double
%                          .psi         backbone dihedral psi double
%              coor        .xyz         coordinates           (N,3) double
%                          .elements    atomic numbers        (N,1) uint8
%                          .indices     atom indices          (N,1) int
%                          .all_indices full index array      (N,5) int16
%
% OUTPUT
% argout       cell array of outputs, see above
% exceptions   cell vector of MException objects if something went wrong, 
%              defaults to one cell holding an empty array
%
% access to a DSSP exectuable is required for attribute dssp

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2020: Gunnar Jeschke

% initialize empty output
argout = {[]};
exceptions = {[]};

% set default arguments
if ~exist('conformer','var') || isempty(conformer)
    conformer = 1;
end

switch attribute
    case 'dssp'
        options.order = conformer;
        fname = sprintf('tmp_%s.pdb',datestr(datetime,'yyyy-mmm-dd.HH.MM.SS'));
        exceptions = put_pdb(entity,fname,options);
        warnings = length(exceptions);
        if isempty(exceptions{1})
            warnings = 0;
        else
            my_exception = exceptions{1};
            if strcmp(my_exception.identifier,'put_pdb:unknown_conformer')
                return
            end
        end
        dssp0 = get_dssp(fname);
        chainfields = fieldnames(entity);
        cpoi = 0;
        for kc = 1:length(chainfields)
            chain = chainfields{kc};
            if isstrprop(chain(1),'upper')
                cpoi = cpoi + 1;
                chainfields{cpoi} = chainfields{kc};
            end
        end
        chainfields = chainfields(1:cpoi);
        % residue pointers for all chains
        pointers = zeros(1,100);
        for k = 1:length(dssp0)
            if isstrprop(dssp0(k).chain,'lower')
                chainfield = strcat(upper(chain),'_');
            else
                chainfield = dssp0(k).chain;
            end
            chain_exists = false;
            for kc = 1:length(chainfields)
                if strcmp(chainfield,chainfields{kc})
                    chain_exists = true;
                end
            end
            if chain_exists
                if ~exist('dssp','var')
                    dssp(1).chain = dssp0(k).chain;
                    cnum = 1;
                    % initialize strings and arrays
                    sequence = pad('?',9999);
                    secondary = sequence;
                    res_numbers = zeros(1,9999);
                    sheets = zeros(9999,2);
                    bp = zeros(9999,2);
                    acc = zeros(1,9999);
                    NHO(9999).hp = [];
                    NHO(9999).en = [];
                    OHN(9999).hp = [];
                    OHN(9999).en = [];
                    tco = zeros(1,9999);
                    kappa = zeros(1,9999);
                    alpha = zeros(1,9999);
                    phi = zeros(1,9999);
                    psi = zeros(1,9999);                 
                    pointers(cnum) = 0;
                else
                    assigned = false;
                    for kc = 1:length(dssp)
                        if strcmp(dssp(kc).chain,dssp0(k).chain)
                            cnum = kc;
                            assigned = true;
                        end
                    end
                    if ~assigned
                        % store information for previous chain
                        dssp(cnum).sequence = sequence(1:pointers(cnum)); %#ok<AGROW>
                        dssp(cnum).secondary = secondary(1:pointers(cnum)); %#ok<AGROW>
                        dssp(cnum).resnum = res_numbers(1:pointers(cnum)); %#ok<AGROW>
                        dssp(cnum).sheets = sheets(1:pointers(cnum),:); %#ok<AGROW>
                        dssp(cnum).bp = bp(1:pointers(cnum),:); %#ok<AGROW>
                        dssp(cnum).acc = acc(1:pointers(cnum)); %#ok<AGROW>
                        dssp(cnum).NHO = NHO(1:pointers(cnum)); %#ok<AGROW>
                        dssp(cnum).OHN = OHN(1:pointers(cnum)); %#ok<AGROW>
                        dssp(cnum).tco = tco(1:pointers(cnum)); %#ok<AGROW>
                        dssp(cnum).kappa = kappa(1:pointers(cnum)); %#ok<AGROW>
                        dssp(cnum).alpha = alpha(1:pointers(cnum)); %#ok<AGROW>
                        dssp(cnum).phi = phi(1:pointers(cnum)); %#ok<AGROW>
                        dssp(cnum).psi = psi(1:pointers(cnum)); %#ok<AGROW>
                        cnum = length(dssp)+1;
                        % initialize strings and arrays
                        sequence = pad('?',9999);
                        secondary = sequence;
                        res_numbers = zeros(1,9999);
                        sheets = zeros(9999,2);
                        bp = zeros(9999,2);
                        acc = zeros(1,9999);
                        NHO(9999).hp = [];
                        NHO(9999).en = [];
                        OHN(9999).hp = [];
                        OHN(9999).en = [];
                        tco = zeros(1,9999);
                        kappa = zeros(1,9999);
                        alpha = zeros(1,9999);
                        phi = zeros(1,9999);
                        psi = zeros(1,9999);
                        pointers(cnum) = 0;
                        dssp(cnum).chain = dssp0(k).chain; %#ok<AGROW>
                    end
                end
                % store information for this residue
                pointers(cnum) = pointers(cnum) + 1;
                poi = pointers(cnum);
                res_numbers(poi) = str2double(dssp0(k).tag);
                sequence(poi) = dssp0(k).slc;
                secondary(poi) = dssp0(k).sec;
                sheets(poi,:) = dssp0(k).sheet;
                bp(poi,:) = dssp0(k).bp;
                acc(poi) = dssp0(k).acc;
                NHO(poi).hp = dssp0(k).NHO.hp; %#ok<AGROW>
                OHN(poi).hp = dssp0(k).OHN.hp; %#ok<AGROW>
                NHO(poi).en = dssp0(k).NHO.energy; %#ok<AGROW>
                OHN(poi).en = dssp0(k).OHN.energy; %#ok<AGROW>
                tco(poi) = dssp0(k).tco;
                kappa(poi) = dssp0(k).kappa;
                alpha(poi) = dssp0(k).alpha;
                phi(poi) = dssp0(k).phi;
                psi(poi) = dssp0(k).psi;
            end % if this is a valid chain
        end % loop over all residues in DSSP structure
        if exist('cnum','var')
            % store information for previous chain
            dssp(cnum).sequence = sequence(1:pointers(cnum)); 
            dssp(cnum).secondary = secondary(1:pointers(cnum)); 
            dssp(cnum).resnum = res_numbers(1:pointers(cnum)); 
            dssp(cnum).sheets = sheets(1:pointers(cnum),:); 
            dssp(cnum).bp = bp(1:pointers(cnum),:); 
            dssp(cnum).acc = acc(1:pointers(cnum)); 
            dssp(cnum).NHO = NHO(1:pointers(cnum)); 
            dssp(cnum).OHN = OHN(1:pointers(cnum)); 
            dssp(cnum).tco = tco(1:pointers(cnum)); 
            dssp(cnum).kappa = kappa(1:pointers(cnum)); 
            dssp(cnum).alpha = alpha(1:pointers(cnum)); 
            dssp(cnum).phi = phi(1:pointers(cnum)); 
            dssp(cnum).psi = psi(1:pointers(cnum)); 
            argout = {dssp};
        else
            warnings = warnings + 1;
            exceptions{warnings} = MException('get_conformer:no_dssp_information', 'No DSSP information found for conformer %i',conformer);
        end
        try
            delete(fname);
        catch exception
            warnings = warnings + 1;
            exceptions{warnings} = exception;
        end
    case 'coor'
        [m,~] = size(entity.index_array);
        all_indices = 1:m;
        coor.indices = all_indices(entity.index_array(:,4) == conformer);
        coor.all_indices = entity.index_array(coor.indices,:);
        coor.xyz = entity.xyz(coor.indices,:);
        if max(coor.indices) > length(entity.elements)
            coor.elements = entity.elements(coor.indices-(conformer-1)*length(entity.elements));
        else
            coor.elements = entity.elements(coor.indices);
        end
        argout = {coor};
    otherwise
        exceptions = {MException('get_conformer:unsupported_attribute', 'Attribute %s not supported',attribute)};
end

