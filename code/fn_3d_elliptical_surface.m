function [vtcs, fcs] = fn_3d_elliptical_surface(centre, axis_dir, maj_dir, maj_rad, min_rad, varargin)
%USAGE
%   [vtcs, fcs] = fn_3d_elliptical_surface(pt1, pt2, rad [, n_ang_divisions])
%AUTHOR
%   Paul Wilcox and Co-pilot (2026)
%SUMMARY
%   Create a 3D elliptical surface described by vertices and triangular faces
%INPUTS
%   centre - [1x3] vector of coordinates of centre
%   axis_dir, maj_dir - [1x3] vector of axes directions
%   maj_rad, min_rad - major and minor radi
%   [n_ang_divisions - number of subdivisions in angle, default = 72]
%OUTPUTS
%   vtcs - [n_vtcs x 3] matrix of vertices
%   fcs - [n_fcs x 3] matrix of vertex indices for each face
%NOTES
%   You can plot the result directly using
%   trisurf(fcs, vtcs(:,1), vtcs(:,2), vtcs(:,3));
%--------------------------------------------------------------------------
if numel(varargin) >= 1
    n_ang_divisions = varargin{1};
else
    n_ang_divisions = 72; % default subdivisions
end

% ---- Input checks ----
validateattributes(centre,   {'numeric'}, {'vector','numel',3,'real','finite'});
validateattributes(axis_dir, {'numeric'}, {'vector','numel',3,'real','finite'});
validateattributes(maj_dir, {'numeric'}, {'vector','numel',3,'real','finite'});
validateattributes(maj_rad,      {'numeric'}, {'scalar','real','finite','positive'});
validateattributes(min_rad,      {'numeric'}, {'scalar','real','finite','positive'});
validateattributes(n_ang_divisions, {'numeric'}, {'scalar','integer','>=',3});

centre   = double(centre(:).');    % ensure 1x3 double
axis_dir = double(axis_dir(:).');
maj_dir = double(maj_dir(:).');

%normalise and get minor axis 
axis_dir = axis_dir / norm(axis_dir);
maj_dir = maj_dir / norm(maj_dir);
min_dir = cross(axis_dir, maj_dir);


%start with a disk
[vtcs, fcs] = fn_3d_disk_surface(centre, axis_dir, maj_rad, n_ang_divisions);

%stretch
for i = 1: size(vtcs, 1)
    r = vtcs(i, :) - centre;
    %project onto major and minor axes
    u = dot(r, maj_dir);
    v = dot(r, min_dir);
    %reduce minor axis dim and transform back
    vtcs(i, :) = u * maj_dir + v * min_dir * (min_rad / maj_rad) + centre;
end

    
end
