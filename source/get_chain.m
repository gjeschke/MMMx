function [argout,exceptions] = get_chain(entity,attribute,address)
%
% GET_CHAIN Retrieves attributes of a chain (or molecule)
%
%   [argout,exceptions] = GET_CHAIN(entity,attribute)
%   Provides attribute values and possibly exceptions for chains selected
%   in an entity
%
%   [argout,exceptions] = GET_CHAIN(entity,attribute,address)
%   Provides attribute values and possibly exceptions for chains selected
%   by address
%
% INPUT
% entity       entity in an MMMx format, must be provided
% address      MMMx address, 'selected' refers to the current selection
%              defaults to 'selected'
% attribute    chain attribute to be retrieved, string, defaults to
%              'info'
%              attribute   output                           Type
%              ------------------------------------------------------------
%              info        .name        chain tag           string
%                          .type        chain type          string
%                                       'peptide','DNA',
%                                       'RNA', 'HET',
%                                       'peptide+', 'DNA+',
%                                       'RNA+', 'water'
%                          .sequence    single-letter       string
%                          .resnums     residue numbers     double vector
%              populations conformer populations            (C,1) double
%
% OUTPUT
% argout       cell array of outputs, see above
% exceptions   cell vector of MException objects if something went wrong, 
%              defaults to one cell holding an empty array
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2020: Gunnar Jeschke

% initialize empty output
argout = cell(1,10000);
exceptions = {[]};

% set default arguments
if ~exist('address','var')
    address = 'selected';
end
if ~exist('attribute','var')
    attribute = 'info';
end

% select the objects by the provided address
entity = select(entity,address,true);

if strcmp(attribute,'info')
    d = load('monomer_attributes.mat');
end

outputs = 0; % counter for the number of outputs
index_vector = zeros(1,5,'uint16');
% scan entity for selected locations
chains = fieldnames(entity);
for kc = 1:length(chains)
    chain = chains{kc};
    if isstrprop(chain(1),'upper') % chain fields start with a capital
        index_vector(1) =  entity.(chain).index;
        if entity.(chain).selected
            outputs = outputs + 1;
            switch attribute
                case 'info'
                    cname = chain;
                    if cname(end) == '_'
                        cname(1) = lower(cname(1));
                        cname = cname(1:end-1);
                    end
                    info.name = cname;
                    peptide = false;
                    DNA = false;
                    RNA = false;
                    HET = false;
                    water = false;
                    residues = fieldnames(entity.(chain));
                    resnums = zeros(length(residues),1);
                    sequence = pad('?',length(residues));
                    nres = 0;
                    for kr = 1:length(residues) % expand over all residues
                        residue = residues{kr};
                        if strcmp(residue(1),'R') % these are residue fields
                            nres = nres + 1;
                            resnums(nres) = str2double(residue(2:end));
                            rtag = entity.(chain).(residue).name;
                            res = tag2id(rtag,upper(d.monomers.aa_tags));
                            if ~isempty(res)
                                peptide = true;
                                res = tag2id(rtag,upper(d.monomers.aa_tags));
                                sequence(nres) = d.monomers.aa_slc(res);
                            elseif contains(upper(d.monomers.nt_tags),rtag)
                                if rtag(1) == 'D'
                                    DNA = true;
                                    rtag = [' ' rtag]; %#ok<AGROW>
                                else
                                    RNA = true;
                                    rtag = ['  ' rtag]; %#ok<AGROW>
                                end
                                res = tag2id(rtag,upper(d.monomers.nt_tags));
                                sequence(nres) = d.monomers.nt_slc(res);
                            elseif contains(upper(d.monomers.nt_tags_CYANA),rtag)
                                if rtag(1) == 'R' || rtag(1) == 'U'
                                    RNA = true;
                                else
                                    DNA = true;
                                end
                                res = tag2id(rtag,upper(d.monomers.nt_tags_CYANA));
                                sequence(nres) = d.monomers.nt_slc(res);
                            elseif strcmp(rtag,'HOH')
                                water = true;
                                sequence(nres) = '!';
                            else
                                HET = true;
                                sequence(nres) = '?';
                            end
                        end
                    end
                    mytype = 'water';
                    if peptide
                        mytype = 'peptide';
                        if HET || DNA || RNA || water
                            mytype = 'petide+';
                        end
                    elseif RNA
                        mytype = 'RNA';
                        if DNA || HET || water
                            mytype = 'RNA+';
                        end
                    elseif DNA
                        mytype = 'DNA';
                        if HET || water
                            mytype = 'DNA+';
                        end
                    elseif HET
                        mytype = 'HET';
                    end
                    info.type = mytype;
                    info.sequence = sequence(1:nres);
                    info.resnums = resnums(1:nres);
                    argout{outputs} = info;
                case 'populations'
                    argout{outputs} = entity.populations;
                otherwise
                    argout = {};
                    exceptions = {MException('get_chain:unsupported_attribute', 'Attribute %s not supported',attribute)};
                    return
            end
        end
    end
end

argout = argout(1:outputs);

