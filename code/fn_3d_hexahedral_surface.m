function [vtcs, fcs] = fn_3d_hexahedral_surface(pts)
%USAGE
%   [vtcs, fcs] = fn_3d_hexahedral_surface(pts)
%AUTHOR
%   Paul Wilcox (2026)
%SUMMARY
%   Create a 3D hexahedral surface described by vertices and triangular 
%   faces based list of 8 corner vertices, ordered cyclically around 
%   'bottom' face and the cyclically around 'top' face in same sense (i.e.
%   like Abaqus node numbering)
%INPUTS
%   pts - [8x3] vector of coordinates corner vertices
%OUTPUTS
%   vtcs - [n_vtcs x 3] matrix of vertices
%   fcs - [n_fcs x 3] matrix of vertex indices for each face
%NOTES
%   You can plot the result directly using
%   trisurf(fcs, vtcs(:,1), vtcs(:,2), vtcs(:,3));
%--------------------------------------------------------------------------

% ---- Input checks ----
validateattributes(pts, {'numeric'}, {'size', [8, 3], 'real','finite'});

pts = double(pts);   % ensure 1x3 double

vtcs = pts;
fcs = fn_polygons_to_triangles([
                1,2,3,4
                1,2,6,5
                2,3,7,6
                3,4,8,7
                4,1,5,8
                5,6,7,8]);

end