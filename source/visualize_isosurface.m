function visualize_isosurface(density_file,property_file,options)
% visualize_isosurface(density_file,property_file,options)
%
% Visualizes volume data as an isosurface, which can be colored by property
% volume data and can be saved to a figure file as an option
% further options allow control of the visualization
% 
% INPUT
% density_file  name of a Matlab .mat file that contains a 3D array 'cube'
%               and corresponding axes 'x', 'y', and 'z'
%               the isosurface is constructed from these volume data
% property_file name of a Matlab .mat file that contains a 3D array 'cube'
%               and corresponding axes 'x', 'y', and 'z', optional
%               these volume data determine coloring, if present
% options       control options, structure with fields
%               .level          determines isovalue, number between 0 and 1,
%                               is the fraction of the volume included by 
%                               the isosurface, defaults to 0.75
%               .camvec         camera position vector, with the origin at 
%                               the center of the plot, can be 'x', '-x',
%                               'y', '-y', 'z', '-z' or a vector [1,3],
%                               defaults to 'x'
%               .colorscheme    flat color of the isosurface or color
%                               scheme if a property file is present, for a
%                               flat color, this can be a SVG color name or
%                               an RGB triple, default is 'gainsboro'
%                               for a color scheme, this can be
%                               'electrostatic', 'cation-pi', or
%                               'hydrophobic', with the Matlab parula
%                               colormap as default
%               .opaqueness     opaqueness of the isosurface, value between
%                               0 and 1, with 1 being dully opaque and 0
%                               fully transparent, defaults to 1
%               .limits         used only with property file, limits of the
%                               color axis, defaults to [minimum, maximum]
%                               of the face color data unless the color
%                               scheme is
%                               electrostatic : -1,1
%                               cation-pi     : -0.2,0.2
%                               hydrophobic   : -10,10
%                               can be a vector [min,max] or a single value
%                               scale, then corresponding to [-scale,scale]
%               .camupvec       upward direction of camera, vector [1,3],
%                               defaults to [0,0,1], must not be parallel
%                               to camvec
%               .background     figure background color, can be an SVG
%                               color name or an RGB triple, not used when
%                               the figure is copied or saved, edit ccopy
%                               options manually, if you want to save with
%                               non-transparent background
%               .figname        file name for saving the figure, figure is
%                               saved automatically only .figname is
%                               present, the extension determines file
%                               format, .png is recommended

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2023: Gunnar Jeschke

% visualize_isosurface(density_file,level,camvec,colorscheme,opaqueness,property_file,limits,camupvec,background,figname)

if ~exist('options','var') || isempty(options)
    options.level = 0.75;
end

if ~isfield(options,'camvec')
    options.camvec = [1,0,0];
elseif ischar(options.camvec)
    switch options.camvec
        case 'x'
            options.camvec = [1,0,0];
        case '-x'
            options.camvec = [-1,0,0];
        case 'y'
            options.camvec = [0,1,0];
        case '-y'
            options.camvec = [0,-1,0];
        case 'z'
            options.camvec = [0,0,1];
        case '-z'
            options.camvec = [0,0,-1];
        otherwise
            options.camvec = [1,0,0];
    end
end

if ~isfield(options,'camupvec')
    options.camupvec = [0,0,1];
end

if ~isfield(options,'colorscheme') || isempty(options.colorscheme)
    if ~exist('property_file','var') || isempty(property_file)
        options.colorscheme = 'gainsboro';
    else
        options.colorscheme = 'electrostatic';
    end
end

if ~isfield(options,'opaqueness')
    options.opaqueness = 1;
end

if isfield(options,'limits') && length(options.limits) == 1 && ~isnan(options.limits)
    options.limits = [-options.limits,options.limits];
end

if ~isfield(options,'background')
    options.background = [0.94,0.94,0.94];
end

if ischar(options.background)
    options.background = get_svg_color(options.background);
end

d = load(density_file);
density = d.cube;

total_density = sum(sum(sum(density)));
isovalue = max(max(max(density)));
while sum(sum(sum(density(density >= isovalue)))) < options.level*total_density
    isovalue = 0.99*isovalue;
end

h = figure; clf; 
h.Color = options.background;
if exist('property_file','var') && ~isempty(property_file)
    p = load(property_file);
    s = isosurface(d.y,d.x,d.z,density,isovalue,p.cube);
    patch(s,'FaceColor','interp','EdgeColor','none','FaceLighting','gouraud','BackFaceLighting','lit','FaceAlpha',options.opaqueness);
    fprintf(1,'Range: %6.3f, %6.3f\n',min(s.facevertexcdata),max(s.facevertexcdata));
    if ~isempty(options.limits) && isnan(options.limits)
        options.limits = [-max(s.facevertexcdata),max(s.facevertexcdata)];
        if -min(s.facevertexcdata) > max(s.facevertexcdata)
            options.limits = [min(s.facevertexcdata),-min(s.facevertexcdata)];
        end
    end
    switch options.colorscheme
        case 'cation-pi'
            colormap(colorscale('gold','blue'));
            if ~isfield(options,'limits') || isempty(options.limits)
                options.limits = [-0.2,0.2];
            end
        case 'electrostatic' 
            colormap(colorscale('red','blue'));
            if ~isfield(options,'limits') || isempty(options.limits)
                options.limits = [-1,1];
            end
        case 'hydrophobic'
            colormap(colorscale('green','darkviolet'));
            if ~isfield(options,'limits') || isempty(options.limits)
                options.limits = [-10,10];
            end
        otherwise
            colormap(parula);
    end
else
    s = isosurface(d.y,d.x,d.z,density,isovalue);
    if ischar(options.colorscheme)
        patchcolor = get_svg_color(options.colorscheme);
    else
        patchcolor = options.colorscheme;
    end
    patch(s,'FaceColor',patchcolor,'EdgeColor','none','FaceLighting','gouraud','BackFaceLighting','lit','FaceAlpha',options.opaqueness);
end

my_axes = gca;
if isfield(options,'limits') && ~isempty(options.limits)
    my_axes.CLim = options.limits;
end
view(my_axes,options.camvec);
my_axes.CameraUpVector = options.camupvec;
material shiny
camlight
axis equal
axis off

if isfield(options,'figname')
    if isempty(strfind(options.figname,'.'))
        options.figname = strcat(options.figname,'.png');
    end
    saveas(h,options.figname);
end