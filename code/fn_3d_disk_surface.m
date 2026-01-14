function [vtcs, fcs] = fn_3d_disk_surface(centre, axis_dir, rad, varargin)
%USAGE
%   [vtcs, fcs] = fn_3d_disk_surface(pt1, pt2, rad [, n_ang_divisions])
%AUTHOR
%   Paul Wilcox and Co-pilot (2026)
%SUMMARY
%   Create a 3D disk surface described by vertices and triangular faces
%INPUTS
%   centre - [1x3] vector of coordinates of centre
%   axis_dir - [1x3] vector of axis direction
%   rad - radius
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
    validateattributes(rad,      {'numeric'}, {'scalar','real','finite','positive'});
    validateattributes(n_ang_divisions, {'numeric'}, {'scalar','integer','>=',3});

    centre   = double(centre(:).');    % ensure 1x3 double
    axis_dir = double(axis_dir(:).');

    % ---- Normalize axis and build orthonormal in-plane basis (u, v) ----
    k = axis_dir / norm(axis_dir);
    assert(all(isfinite(k)) && any(k ~= 0), 'axis_dir must be a non-zero vector.');

    % choose a helper vector not parallel to k
    if abs(k(1)) < 0.9
        a = [1 0 0];
    else
        a = [0 1 0];
    end
    u = cross(k, a);  u = u / norm(u);     % first in-plane unit vector
    v = cross(k, u);                        % second in-plane unit vector (unit, orthogonal)

    % ---- Angular samples (rim) ----
    n = n_ang_divisions;
    theta = (0:n-1).' * (2*pi/n);          % column vector of angles

    % direction vectors around the circle (implicit expansion, R2016b+)
    C = cos(theta) * u + sin(theta) * v;   % Nx3

    % ---- Rim vertices and center ----
    rim   = centre + rad * C;              % Nx3
    vtcs  = [rim; centre];                 % (n+1) x 3
    c_idx = size(vtcs,1);                  % center vertex index

    % ---- Triangular facets (fan) with +k outward normal ----
    i1 = (1:n).';
    i2 = mod(i1, n) + 1;                   % wrap-around next index

    % CCW when viewed along +k: [rim(i), rim(i_next), center]
    fcs = [i1, i2, repmat(c_idx, n, 1)];   % n x 3
end
