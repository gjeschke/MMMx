function [coor,indices,elements,exceptions] = get_coor(entity,address,heavy,paradigm)
%
% GET_COOR Retrieves atom location coordinates for selected objects 
%
%   [coor,indices,elements,exceptions] = GET_COOR(entity)
%   Returns a Cartesian coordinate array, indices into the corresponding
%   entity array, elements, and possibly exceptions for selected objects in
%   an entity
%
%   [coor,indices,elements,exceptions] = GET_COOR(entity,address)
%   Returns a Cartesian coordinate array, indices into the corresponding
%   entity array, elements, and possibly exceptions for objects in an
%   entity slected by address
%
%   [coor,indices,elements,exceptions] = GET_COOR(entity,address,heavy)
%   Neglects hydrogen atoms, if heavy is true, heavy defaults to false
%
%   [coor,indices,elements,exceptions] = GET_COOR(entity,address,heavy,paradigm)
%   Returns only the first location(rotamer) in the first conformer for
%   each atom
%
% INPUT
% entity       entity in an MMMx format, must be provided
% address      MMMx address, 'selected' refers to the current selection
%              defaults to 'selected'
% heavy        flag, if true, hydrogen atoms are excluded, defaults to
%              false
% paradigm     flag, if true, only first location of each atom is returned,
%              defaults to false
%
% OUTPUT
% coor         Cartesian coordinates, (N,3) double array
% indices      indices into entity.xyz, entity.elements; (N,1) int array 
% exceptions   cell vector of MException objects if selection by address
%              went wrong, defaults to one cell holding an empty array
% elements     atomic numbers; (N,1) int8 array
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2020: Gunnar Jeschke

% initialize empty output
coor = [];
indices = [];
elements = [];

% set default arguments
if ~exist('address','var') || isempty(address)
    address = 'selected';
end

if ~exist('heavy','var') || isempty(heavy)
    heavy = false;
end

if ~exist('paradigm','var') || isempty(paradigm)
    paradigm = false;
end


% select the objects by the provided address
[entity,exceptions] = select(entity,address,true);

% return and report back exceptions, if selection went wrong

if ~isempty(exceptions) && ~isempty(exceptions{1})
    return
end

% get atom indices
if paradigm % only the first location/rotamer of each atom
    [indices,complete] = get_selection(entity);
    indices = indices((complete(:,4) == 1) & (complete(:,5) == 1));
else
    indices = get_selection(entity);
end

coor = entity.xyz(indices,:);

% retrieve elements, if requested
if nargout >= 3 || heavy
    elements = entity.elements(indices);
end

if heavy
    coor = coor(elements > 1,:);
    indices = indices(elements > 1);
    elements = elements(elements > 1);
end

