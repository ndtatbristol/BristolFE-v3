function [vtcs, fcs] = fn_3d_spherical_surface(cent, rad, varargin)
%USAGE
%   [vtcs, fcs] = fn_3d_spherical_surface(cent, rad [, n_sub_divisions])
%AUTHOR
%   Paul Wilcox + Co-pilot (2025)
%SUMMARY
%   Create a 3D sphere described by vertices and faces
%INPUTS
%   cent - [1x3] vector of coordinates of centre
%   rad - radius
%   [n_sub_divisions - number of subdivisions of icosohedredon to use for output,
%   default = 2, which gives 320 faces
%OUTPUTS
%   vtcs - [n_vtcs x 3] matrix of vertices
%   fcs - [n_fcs x 3] matrix of vertex indices for each face
%NOTES
%   You can plot the result directly using
%   trisurf(fcs, vtcs(:,1), vtcs(:,2), vtcs(:,3));
%--------------------------------------------------------------------------
if numel(varargin) >= 1
    n_sub_divisions = varargin{1};
else
    n_sub_divisions = 2; % default subdivisions
end

% Start from an icosahedron
[vtcs, fcs] = icosahedron();

% Subdivide faces
for i = 1:n_sub_divisions
    [vtcs, fcs] = subdivide(vtcs, fcs);
end

% Normalize vertices to unit sphere, scale and shift
vtcs = vtcs ./ vecnorm(vtcs, 2, 2) * rad + cent;
end

function [V, F] = icosahedron()
phi = (1 + sqrt(5)) / 2; % golden ratio
V = [
    -1,  phi, 0;
    1,  phi, 0;
    -1, -phi, 0;
    1, -phi, 0;
    0, -1,  phi;
    0,  1,  phi;
    0, -1, -phi;
    0,  1, -phi;
    phi, 0, -1;
    phi, 0,  1;
    -phi, 0, -1;
    -phi, 0,  1
    ];
V = V ./ norm(V(1,:)); % normalize initial vertices
F = [
    1,12,6; 1,6,2; 1,2,8; 1,8,11; 1,11,12;
    2,6,10; 6,12,5; 12,11,3; 11,8,7; 8,2,9;
    4,10,5; 4,5,3; 4,3,7; 4,7,9; 4,9,10;
    5,10,6; 3,5,12; 7,3,11; 9,7,8; 10,9,2
    ];
end

function [Vnew, Fnew] = subdivide(V, F)
edgeMap = containers.Map('KeyType','char','ValueType','int32');
Vnew = V;
Fnew = zeros(size(F,1)*4, 3);
nextIdx = size(V,1) + 1;

for i = 1:size(F,1)
    tri = F(i,:);
    mids = zeros(1,3);
    for e = 1:3
        a = tri(e);
        b = tri(mod(e,3)+1);
        key = sprintf('%d-%d', min(a,b), max(a,b));
        if isKey(edgeMap, key)
            mids(e) = edgeMap(key);
        else
            mid = (V(a,:) + V(b,:)) / 2;
            Vnew(nextIdx,:) = mid;
            edgeMap(key) = nextIdx;
            mids(e) = nextIdx;
            nextIdx = nextIdx + 1;
        end
    end
    % Create 4 new triangles
    Fnew(4*(i-1)+1,:) = [tri(1), mids(1), mids(3)];
    Fnew(4*(i-1)+2,:) = [tri(2), mids(2), mids(1)];
    Fnew(4*(i-1)+3,:) = [tri(3), mids(3), mids(2)];
    Fnew(4*(i-1)+4,:) = [mids(1), mids(2), mids(3)];
end
end
