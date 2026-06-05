function [vtcs, fcs] = fn_3d_rectalinear_surface(pt1, pt2)
%USAGE
%   [vtcs, fcs] = fn_3d_rectalinear_surface(pt1, pt2)
%AUTHOR
%   Paul Wilcox (2026)
%SUMMARY
%   Create a 3D rectalinear surface described by vertices and triangular 
%   faces based on two opposing corner positions
%INPUTS
%   pt1, pt2 - [1x3] vector of coordinates of opposing corners
%OUTPUTS
%   vtcs - [n_vtcs x 3] matrix of vertices
%   fcs - [n_fcs x 3] matrix of vertex indices for each face
%NOTES
%   You can plot the result directly using
%   trisurf(fcs, vtcs(:,1), vtcs(:,2), vtcs(:,3));
%--------------------------------------------------------------------------

% ---- Input checks ----
validateattributes(pt1, {'numeric'}, {'vector','numel',3,'real','finite'});
validateattributes(pt2, {'numeric'}, {'vector','numel',3,'real','finite'});

pt1 = double(pt1(:).');   % ensure 1x3 double
pt2 = double(pt2(:).');

% vtcs = [
%     pt1(1), pt1(2), pt1(3)
%     pt2(1), pt1(2), pt1(3)
%     pt1(1), pt2(2), pt1(3)
%     pt2(1), pt2(2), pt1(3)
%     pt1(1), pt1(2), pt2(3)
%     pt2(1), pt1(2), pt2(3)
%     pt1(1), pt2(2), pt2(3)
%     pt2(1), pt2(2), pt2(3)];

% fcs = [
%     1,2,3 %bottom1
%     2,4,3 %bottom2
%     1,2,5 %front1
%     2,6,5 %front2
%     2,4,6 %right1
%     4,8,6 %right2
%     3,1,7 %left1
%     1,5,7 %left2
%     4,3,7 %back1
%     8,4,7 %back2
%     5,6,7 %top1
%     6,8,7 %top2
%     ];

vtcs = [
    pt1(1), pt1(2), pt1(3)
    pt2(1), pt1(2), pt1(3)
    pt2(1), pt2(2), pt1(3)
    pt1(1), pt2(2), pt1(3)
    pt1(1), pt1(2), pt2(3)
    pt2(1), pt1(2), pt2(3)
    pt2(1), pt2(2), pt2(3)
    pt1(1), pt2(2), pt2(3)];

fcs = fn_polygons_to_triangles([
                1,2,3,4
                1,2,6,5
                2,3,7,6
                3,4,8,7
                4,1,5,8
                5,6,7,8]);

end