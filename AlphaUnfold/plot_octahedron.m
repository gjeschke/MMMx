function obj = plot_octahedron(xyz, a, b, c, color)
% PLOT_OCTAHEDRON Plots an octahedron with opposite vertices along x, y, z axes
%   plot_octahedron(a, b, c, color) plots an octahedron where:
%       xyz = coordinate of the center
%       a = distance between opposite vertices along x-axis
%       b = distance between opposite vertices along y-axis
%       c = distance between opposite vertices along z-axis
%       color = 1x3 RGB vector specifying face color (e.g., [0.5 0.5 1])

% Define the vertices
% Vertices are arranged with opposite pairs along each axis
vertices = [
    a/2,  0,    0;    % Vertex 1: +x
   -a/2,  0,    0;    % Vertex 2: -x
    0,   b/2,   0;    % Vertex 3: +y
    0,  -b/2,   0;    % Vertex 4: -y
    0,    0,   c/2;    % Vertex 5: +z
    0,    0,  -c/2     % Vertex 6: -z
];
vertices = vertices + repmat(xyz,6,1);

% Define the faces (each row is a triangular face with 3 vertex indices)
% The octahedron has 8 triangular faces
faces = [
    1, 3, 5;    % Face 1: +x, +y, +z
    1, 5, 4;    % Face 2: +x, +z, -y
    1, 4, 6;    % Face 3: +x, -y, -z
    1, 6, 3;    % Face 4: +x, -z, +y
    2, 5, 3;    % Face 5: -x, +z, +y
    2, 3, 6;    % Face 6: -x, +y, -z
    2, 6, 4;    % Face 7: -x, -z, -y
    2, 4, 5     % Face 8: -x, -y, +z
];

% Create the patch object
obj = patch('Faces', faces, 'Vertices', vertices);
set(obj, 'FaceColor', color, 'EdgeColor', 'none',...
    'FaceAlpha',1,'FaceLighting','gouraud','Clipping','off');
set(obj, 'CDataMapping','direct','AlphaDataMapping','none');


% Set axis properties for better visualization
% axis equal;
% grid on;
% xlabel('X');
% ylabel('Y');
% zlabel('Z');
% title(sprintf('Octahedron: a=%.1f, b=%.1f, c=%.1f', a, b, c));

% Set view for 3D visualization
% view(3);
% rotate3d on;

% Adjust axes limits to fit the octahedron
% margin = 0.1;
% xlim([-a/2*(1+margin), a/2*(1+margin)]);
% ylim([-b/2*(1+margin), b/2*(1+margin)]);
% zlim([-c/2*(1+margin), c/2*(1+margin)]);

end