function d = fn_quick_dist_to_point_bdry(pts, bdry_nds)
%SUMMARY
%   High speed version of fn_{2d/3d}_signed_dist_to_bdry that calculates
%   minimum distance between each point in pts and all bdry_nds. Distance
%   is not signed, and function is primarily for use when calculating
%   absorbing level when boundary is defined by conneced node on FE mesh
%   (i.e. very dense, so loss of accuracy by not considering vertices and
%   facets separately is minimal.

d = ones(size(pts, 1), 1) * inf;

for i = 1:size(bdry_nds, 1)
    r = sum((bdry_nds(i, :) - pts) .^ 2, 2);
    d = min(r, d);
end

d = sqrt(d);
end