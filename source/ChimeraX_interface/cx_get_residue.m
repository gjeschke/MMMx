function [attribute,exceptions] = cx_get_residue(id,address,whichone)
%
% CX_GET_RESIDUE Set an attribute of an atom in ChimeraX
%
%   attribute = CX_GET_RESIDUE(id,address,whichone,attribute)
%   Gets value of an atom atrribute in a ChimeraX entity
%
%   [attribute,exceptions] = CX_GET_RESIDUE(id,address,whichone,attribute)
%
% INPUT
% id        ChimeraX internal identifier for entity
% address   MMMx address for the residue
% whichone  valid attribute name
%           center        (1,3) double center coordinate
%           chi1          double, side chain dihedral 1
%           chi2          double, side chain dihedral 2
%           chi3          double, side chain dihedral 3
%           chi4          double, side chain dihedral 4
%           color         (1,3) int RGB 0... 255
%           hide          flag
%           name          string, residue name (three-letter code)
%           omega         double, backbone angle (cis/trans)
%           phi           double, backbone angle
%           polymer_type  int, 0 none, 1 peptide, 2 nucleic acid
%           psi           double, backbone angle
%           secondary     char, 'H' helix, 'E' strand, '' none 
%           selected      Boolean flag
%           tlc           string, three-letter code
%           transparency  double
%
% OUTPUT
% attribute  attribute value, empty if the attribute is not supported or
%            cannot be retrieved
% exceptions cell vector of MException objects if something went wrong, 
%            defaults to one cell holding an empty array
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2020: Gunnar Jeschke

% ChimeraX requires plural for levels, but MMMx accepts also singular

% initialize empty exception cell vector

% allow for capitalization in attribute name
whichone = lower(whichone);

[spec,exceptions] = cx_from_mmmx_address(address,id);
if ~isempty(exceptions{1})
    return
end

% special treatment for color transparency
if strcmp(whichone,'transparency') || strcmp(whichone,'color')
    attribute = cx_get_attribute(spec,'residue','ribbon_color');
elseif strcmp(whichone,'secondary')
    attribute = '';
    attribute1 = cx_get_attribute(spec,'residue','is_helix');
    if contains(attribute1,'True')
        attribute = 'H';
    end
    attribute2 = cx_get_attribute(spec,'residue','is_strand');
    if contains(attribute2,'True')
        attribute = 'E';
    end
elseif strcmp(whichone,'hide')
    attribute = cx_get_attribute(spec,'residue','ribbon_display');
elseif strcmp(whichone,'tlc')
    attribute = cx_get_attribute(spec,'residue','name');
else
    attribute = cx_get_attribute(spec,'residue',whichone);
end


attribute = splitlines(attribute);
attribute = attribute{1};

switch whichone
    case 'center'
        center_str = split(attribute);
        attribute = str2num(center_str{1}); %#ok<ST2NM>
    case {'chi1','chi2','chi3','chi4'}
        chi_str = split(attribute);
        attribute = str2double(chi_str{1});
    case 'color'
        col_str = split(attribute);
        attribute = str2num(col_str{1}); %#ok<ST2NM>
        attribute = attribute(1:3);
    case 'secondary'
    case 'tlc'
        name_str = split(attribute);
        attribute = upper(name_str{1});
    case 'name'
        name_str = split(attribute);
        attribute = upper(name_str{1});
    case {'omega','phi','psi'}
        ang_str = split(attribute);
        attribute = str2double(ang_str{1});
    case 'polymer_type'
        poly_str = split(attribute);
        attribute = str2double(poly_str{1});
    case 'hide'
        disp_str = split(attribute);
        attribute = disp_str{1};
        switch attribute
            case 'False'
                attribute = 1;
            case 'True'
                attribute = 0;
        end
    case 'selected'
        disp_str = split(attribute);
        attribute = disp_str{1};
        switch attribute
            case 'False'
                attribute = 0;
            case 'True'
                attribute = 1;
        end
    case 'transparency'
        col_str = split(attribute);
        attribute = str2num(col_str{1}); %#ok<ST2NM>
        attribute = 100 - 100*attribute(4)/255;
    otherwise
        exceptions{1} = MException('cx_get_residue:no_such_attribute', 'Attribute %s not supported for a residue', whichone);
end