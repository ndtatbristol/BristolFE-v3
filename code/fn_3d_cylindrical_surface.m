function [vtcs, fcs] = fn_3d_cylindrical_surface(pt1, pt2, rad, varargin)
%USAGE
%   [vtcs, fcs] = fn_3d_cylindrical_surface(pt1, pt2, rad [, n_ang_divisions])
%AUTHOR
%   Paul Wilcox and Co-pilot (2026)
%SUMMARY
%   Create a 3D cylindrical surface described by vertices and triangular faces
%INPUTS
%   pt1, pt2 - [1x3] vector of coordinates of end points of axis
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
validateattributes(pt1, {'numeric'}, {'vector','numel',3,'real','finite'});
validateattributes(pt2, {'numeric'}, {'vector','numel',3,'real','finite'});
validateattributes(rad, {'numeric'}, {'scalar','real','finite','positive'});
validateattributes(n_ang_divisions, {'numeric'}, {'scalar','integer','>=',3});

pt1 = double(pt1(:).');   % ensure 1x3 double
pt2 = double(pt2(:).');
axis_vec = pt2 - pt1;
h = norm(axis_vec);
assert(h > 0, 'pt1 and pt2 must be distinct points.');

% ---- Build an orthonormal frame (u, v, k) around the cylinder axis ----
k = axis_vec / h;  % unit axis direction

% pick a helper vector not parallel to k
if abs(k(1)) < 0.9
    a = [1 0 0];
else
    a = [0 1 0];
end
u = cross(k, a);  u = u / norm(u);
v = cross(k, u);  % orthonormal with u and k

% ---- Angular samples ----
n = n_ang_divisions;
theta = (0:n-1).' * (2*pi/n);  % column vector

% circle frame directions per angle (Nx3)
C = cos(theta) * u + sin(theta) * v;  % implicit expansion (R2016b+)

% ---- Two rings of vertices ----
ring1 = pt1 + rad * C;   % at pt1 plane
ring2 = pt2 + rad * C;   % at pt2 plane

vtcs = [ring1; ring2];   % (2n) x 3

% ---- Add center vertices for the two caps ----
idx_c1 = size(vtcs,1) + 1;  % center at pt1
idx_c2 = idx_c1 + 1;        % center at pt2
vtcs   = [vtcs; pt1; pt2];  % (2n+2) x 3

% ---- Side faces as triangles ----
% Indices around the ring
i1 = (1:n).';
i2 = mod(i1, n) + 1;               % wrap-around to next index

% Each quad [i1, i2, i2+n, i1+n] splits into two triangles:
% Triangle 1: [i1, i2, i2+n]
% Triangle 2: [i1, i2+n, i1+n]
side_tri_1 = [i1, i2, i2 + n];     % n x 3
side_tri_2 = [i1, i2 + n, i1 + n]; % n x 3

% ---- End caps as triangle fans with outward normals ----
% Orientation notes:
% - Ring order (increasing theta) is CCW when viewed along +k.
% - Bottom cap (at pt1) outward normal is -k -> use clockwise order when viewed along +k.
% - Top cap (at pt2) outward normal is +k -> use CCW order when viewed along +k.
cap1_tris = [i2, i1, repmat(idx_c1, n, 1)];          % bottom: [i2, i1, c1]
cap2_tris = [i1 + n, i2 + n, repmat(idx_c2, n, 1)];  % top:    [i1+n, i2+n, c2]

% ---- Concatenate all triangles into a single numeric face matrix ----
fcs = [side_tri_1; side_tri_2; cap1_tris; cap2_tris];  % (4n) x 3
end
