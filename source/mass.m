function mass = mass(elements)
%
% MASS Computes mass from an element vector 
%
%   mass = MASS(elements)
%   Adds the mass of all atoms acording to their atomic number in input
%   vector elements, assuming natural isotope abundance
%
%   requires access to defs/element_attributes, returns empty output if
%   this is not accessible
%
% INPUT
% elements     vector of atomic numbers
%
% OUTPUT
% mass         total mass
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2020: Gunnar Jeschke

% initialize empty output
mass = [];

try
    d = load('element_attributes.mat');
catch exception %#ok<NASGU>
    return;
end

masses = d.pse.mass;
masses = masses(elements);
mass = sum(masses);


