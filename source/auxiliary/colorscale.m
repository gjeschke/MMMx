function map = colorscale(negative,positive,grades)
% map = colorscale(negative,positive,grades)
%
% generates a color map from color negative via white to color positive
% with a certain number of grades from one of the extreme colors to white
% useful for electrostatic, cation-pi, and hydrophobic plots
%
% INPUT
%
% negative  either RGB triplet or SVG color name for negative color
% positive  either RGB triplet or SVG color name for positive color
% grades    number of grades for one half, defaults to 25
%
% OUTPUT
%
% map   colormap of dimension [2*grades+1,3]
%
% G. Jeschke, 14.04.2023

if ~exist('grades','var') || isempty(grades)
    grades = 25;
end

if ischar(negative)
    negative = get_svg_color(negative);
end

if ischar(positive)
    positive = get_svg_color(positive);
end

map = ones(2*grades+1,3);
for k = 1:grades+1
    colored = (grades-k+1)/grades;
    map(k,:) = negative*colored + [1,1,1]*(1-colored);
    map(2*grades+1-k+1,:) = positive*colored + [1,1,1]*(1-colored);
end
