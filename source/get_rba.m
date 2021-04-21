function [entity,exceptions] = get_rba(entity,number)
%
% GET_RBA Retrieves a single rigid-body arrangement from an MMMx:RigiFlex
%         model
%
%   entity = GET_RBA(entity,number)
%   Returns a single-conformer entity corresponding to one rigid-body
%   arrangment in a Rigi-generated rigid-nody arrangement ensemble
%
% INPUT
% entity       entity in an MMMx:RigiFlex format, must be provided
% number       number of the rigid-body arrangment to be extracted,
%              defaults to 1
%
% OUTPUT
% entity       output entity, information on other rigid-body arrangements
%              is removed
% exceptions   cell vector of MException objects if something went wrong, 
%              defaults to one cell holding an empty array
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2021: Gunnar Jeschke

% initialize empty output
exceptions = {[]};

% set default arguments
if ~exist('number','var') || isempty(number)
    number = 1;
end

if ~isfield(entity,'rba')
    entity = [];
    exceptions = {MException('get_rba:entity_not_RigiFlex', 'Input entity does not specify rigid-body arrangments')};
    return
end

[nrba,~] = size(entity.rba);
if number < 1 || number > nrba
    entity = [];
    exceptions = {MException('get_rba:requested_rba_missing', 'The requested rigid-body arrangement %i does not exist',number)};
    return
end

nrb = length(entity.rigidbodies); % number of rigid bodies
transvecs = entity.rba(number,:); % transformation vectors

% transform coordinates and remove rba assignment from chains
for kb = 1:nrb
    indices = entity.rigidbodies{kb};
    coor = entity.xyz(indices,:);
    bas = 6*(kb-1);
    trans = transvecs(bas+1:bas+3);
    euler = transvecs(bas+4:bas+6);
    transmat = transrot2affine(trans,euler);
    coor = affine_coor_set(coor,transmat);
    entity.xyz(indices,:) = coor;
    chains = entity.rba_chains{kb};
    for kc = 1:length(chains)
        entity.(chains(kc)) = rmfield(entity.(chains(kc)),'rigidbody');
    end
end

% remove fields that specify all rigid-body arrangments
entity = rmfield(entity,'rba');
entity = rmfield(entity,'rba_populations');
entity = rmfield(entity,'rigidbodies');
entity = rmfield(entity,'rba_chains');

