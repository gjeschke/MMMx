function [entity,exceptions] = set_coor(entity,coor,indices)
%
% SET_COOR Assigns (modified) coordinates to selected atom locations 
%
%   [entity,exceptions] = SET_COOR(entity,coor)
%   Assigns the input coordinates 'coor' to the currently selected objects
%   in entity
%
%   [entity,exceptions] = SET_COOR(entity,coor,indices)
%   Assigns the input coordinates 'coor' to atom locations specified by
%   'indices'
%
%   [entity,exceptions] = SET_COOR(entity,coor,address)
%   Assigns the input coordinates 'coor' to atom locations specified by
%   'address'
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
% entity       the input entity with modified atom coordinates
% exceptions   cell vector of MException objects if selection by address
%              went wrong, defaults to one cell holding an empty array
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2020: Gunnar Jeschke

% initialize empty output
exceptions = {[]};

% set default arguments
if ~exist('indices','var') || isempty(indices)
    address = 'selected';
    indices = [];
end
% if the third argument is a string, interpret it as an address
if ischar(indices)
    address = indices;
    indices = [];
end

% retrieve indices, if not provided
if isempty(indices)
    % select the objects by the provided address
    [entity,exceptions] = select(entity,address,true);
    % return and report back exceptions, if selection went wrong
    if ~isempty(exceptions)
        return
    end
    % get atom indices
    indices = get_selection(entity);
end

entity.xyz(indices,:) = coor;
