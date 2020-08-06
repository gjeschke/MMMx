function [entity,exceptions] = set_residue(entity,attribute,argin,address)
%
% SET_RESIDUE Assigns attributes to residues
%
%   [entity,exceptions] = SET_RESIDUE(entity,attribute,argin)
%   Modifies an attribute of selected residues in entity to the values 
%   provided in argin
%
%   [argout,exceptions] = SET_RESIDUE(entity,attribute,address)
%   Modifies an attribute of residues in entity to the values specified by 
%   address provided in argin
%
% INPUT
% entity       entity in an MMMx format, must be provided
% attribute    residue attribute to be retrieved, string, defaults to
%              'info'
%              attribute   output                           Type
%              ------------------------------------------------------------
%              dssp        DSSP secondary structure         char
%              number      residue number, be supercareful! int
%              populations rotamer populations              (R,1) double
%              sheet       DSSP sheet information           (1,2) double
%              tlc         three-letter code, be careful!   string
% argin        cell array of inputs, see table above
% address      MMMx address, 'selected' refers to the current selection
%              defaults to 'selected'
%
% OUTPUT
% entity       the input entity with  modified attributes
% exceptions   cell vector of MException objects if something went wrong, 
%              defaults to one cell holding an empty array
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2020: Gunnar Jeschke

% initialize empty output
exceptions = {[]};
warnings = 0;

% set default arguments
if ~exist('address','var')
    address = 'selected';
end
if ~exist('attribute','var') || isempty(attribute)
    exceptions = {MException('set_residue:missing_attribute', 'Attribute must be provided')};
end
if ~exist('argin','var')
    exceptions = {MException('set_residue:missing_input', 'Input arguments must be provided')};
end

% select the objects by the provided address
entity = select(entity,address,true);

inputs = 0; % counter for the number of outputs
index_vector = zeros(1,5,'uint16');
% scan entity for selected locations
chains = fieldnames(entity);
for kc = 1:length(chains)
    chain = chains{kc};
    if isstrprop(chain(1),'upper') % chain fields start with a capital
        index_vector(1) =  entity.(chain).index;
        residues = fieldnames(entity.(chain));
        for kr = 1:length(residues) % expand over all residues
            residue = residues{kr};
            if strcmp(residue(1),'R') % these are residue fields
                index_vector(2) =  entity.(chain).(residue).index;
                if entity.(chain).(residue).selected % report back
                    inputs = inputs + 1;
                    switch attribute
                        case 'dssp'
                            entity.(chain).(residue).dssp = argin{inputs};
                        case 'number'
                            resfield = sprintf('R%i',argin{inputs});
                            % check whether this residue number does 
                            % already exist
                            claimed = false;
                            for krtest = 1:length(residues)
                                if strcmp(residues{krtest},resfield)
                                    claimed = true;
                                end
                            end
                            % if so, raise exception and return empty
                            % entity
                            if claimed
                                exceptions = {MException('set_residue:residue_number_exists',...
                                    'Residue number %i does already exist',argin{inputs})};
                                entity = [];
                                return
                            end
                            data = entity.(chain).(residue);
                            modified_chain = rmfield(entity.(chain),residue);
                            modified_chain.(resfield) = data;
                            entity.(chain) = modified_chain;
                        case 'populations'
                            entity.(chain).(residue).populations = argin{inputs};
                        case 'sheet'
                            entity.(chain).(residue).sheet = argin{inputs};
                        case 'tlc'
                            % capitalize
                            tlcinput = upper(argin{inputs});
                            % initialize three-letter code 
                            tlc = char(32*ones(1,3));
                            % remove characters that are not alphanumeric
                            tlcpoi = 0;
                            for kchar = 1:length(tlcinput)
                                if isstrprop(tlcinput(kchar),'alphanum')
                                    tlcpoi = tlcpoi + 1;
                                    tlc(tlcpoi) = tlcinput(kchar);
                                end
                            end
                            % truncate to three characters if necessary
                            if length(tlc) > 3
                                tlc = tlc(1:3);
                            end
                            tlc = strtrim(tlc);
                            % return warning, if three-letter code was
                            % modified
                            if ~strcmp(argin{inputs},tlc)
                                warnings = warnings + 1;
                                exceptions{warnings} = MException('set_residue:three_letter_code_modified',...
                                    'Requested three-letter code %s was modified to %s',argin{inputs},tlc);
                            end
                            % if empty string, raise exception and return
                            % empty entity
                            if isempty(tlc)
                                exceptions = {MException('set_residue:invalid_three_letter_code',...
                                    'Requested three-letter code %s is invalid',tlcinput)};                                
                                entity = [];
                                return
                            end
                            entity.(chain).(residue).name = tlc;
                        otherwise
                            exceptions = {MException('set_residue:unsupported_attribute', 'Attribute %s not supported',attribute)};
                            return
                    end
                end
            end
        end
    end
end

