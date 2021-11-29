function rgb = get_svg_color(name)
% function rgb = = GET_SVG_COLOR(name)
%
% Returns the RGB triple (range 0...1) for an SVG defined color
%
% Input:
%
% name      one of the defined SVG color names, case insensitive, if name
%           is 'help', the output argument is the whole svg color structure
%
% Output:
%
% rgb       rgb triple in range [0,1] for each color, empty if the color
%           name is undefined

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2021: Gunnar Jeschke


colors = load('color_defs');

name = lower(name);
if strcmp(name,'help')
    rgb = colors.svg_colors;
    return
end
    
if ~isfield(colors.svg_colors,name)
    rgb = [];
    return
end

rgb = colors.svg_colors.(name);