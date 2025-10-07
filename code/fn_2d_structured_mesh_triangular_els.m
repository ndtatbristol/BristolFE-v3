function mod = fn_2d_structured_mesh_triangular_els(bdry_pts, el_size, varargin)
%USAGE
%   mod = fn_2d_structured_mesh_triangular_els(bdry_pts, el_size [, force_equilateral_els])
%AUTHOR
%   Paul Wilcox (2025)
%SUMMARY
%   Utility function for generating a isometric structured mesh of triangular
%   elements, that fills the region specified by bdry_nds.
%INPUTS
%   bdry_pts - n_bdry x 2 matrix of coordinates of the n_bdry points that 
%       will define boundary of mesh.
%   el_size - element size
%   [force_equilateral_els = 0] - force elements to be equilaterial 
%   triangles. If bdry elements define a rectangle and 
%   force_equilateral_els = 0, element edges at top and bottom will be
%   exactly on boundary but elements won't be equilateral; if 
%   force_equilateral_els = 1 then element edges may not be exactly on
%   boundary at top and bottom but they will be equilateral.
%OUTPUT
%   mod - structured variable containing fields:
%       .nds - n_nds x 2 matrix of coordinates of each of n_nds nodes
%       .els - n_els x 3 matrix of node indices for each of n_els elements
%       .el_mat_i - n_els x 1 matrix of ones as a placeholder for element
%       material indices assigned elsewhere if more than one type of
%       material is used in model
%       .el_abs_i - n_els x 1 matrix of zeros as a placeholder for element
%       relative absorption
%       .el_typ_i - n_els x 1 matrix of zeros as a placeholder for element
%       type indices
%--------------------------------------------------------------------------
%Figure out a bounding rectangle for whole shape
sin60 = sind(60);

if numel(varargin) > 0
    force_equilateral_els = varargin{1};
else
    force_equilateral_els = 0;
end

if force_equilateral_els
    %Figure out a bounding rectangle for whole shape that is a bit oversize
    %(excess elements will be trimmed later)
    crnr_pts = [min(bdry_pts) - 2 * el_size; max(bdry_pts) + 2 * el_size];
    x = crnr_pts(1, 1):el_size:crnr_pts(2, 1);
    y = crnr_pts(1, 2):el_size * sin60:crnr_pts(2, 2);
else
    %In this case corner points are hard on bdry
    crnr_pts = [min(bdry_pts); max(bdry_pts)];
    %Work out how many nodes are needed in x and y
    nodes_in_x_direction = ceil((crnr_pts(2, 1) - crnr_pts(1, 1)) / el_size) + 1;
    nodes_in_y_direction = ceil((crnr_pts(2, 2) - crnr_pts(1, 2)) / el_size / sin60) + 1;
    x = linspace(crnr_pts(1, 1), crnr_pts(2, 1), nodes_in_x_direction);
    y = linspace(crnr_pts(1, 2), crnr_pts(2, 2), nodes_in_y_direction);
end

%Create cartesian grid of nodal coordinates describing squares
[node_x_positions, node_y_positions] = meshgrid(x, y);

%then shuffle rows of x positions back/forward by half an element to get a
%grid of triangles
el_size_x = x(2) - x(1);
node_x_positions(1:2:end, :) = node_x_positions(1:2:end, :) + el_size_x / 4;
node_x_positions(2:2:end, :) = node_x_positions(2:2:end, :) - el_size_x / 4;

%Work out node numbers associated with each element (a bit fiddly as you can see) 
node_numbers = reshape([1:numel(node_x_positions)], numel(y), numel(x));

element_node1a = node_numbers(1:2:end-1, 1:end-1);
element_node2a = node_numbers(2:2:end, 2:end);
element_node3a = node_numbers(2:2:end, 1:end-1);

element_node1b = node_numbers(1:2:end-1, 1:end-1);
element_node2b = node_numbers(1:2:end-1, 2:end);
element_node3b = node_numbers(2:2:end, 2:end);

element_node1c = node_numbers(2:2:end-1, 2:end);
element_node2c = node_numbers(3:2:end, 2:end);
element_node3c = node_numbers(3:2:end, 1:end-1);

element_node1d = node_numbers(2:2:end-1, 1:end-1);
element_node2d = node_numbers(2:2:end-1, 2:end);
element_node3d = node_numbers(3:2:end, 1:end-1);

%Final m x 2 matrix of x and y coordinates for each node
mod.nds = [node_x_positions(:), node_y_positions(:)];

%Final n x 3 matrix of 3 node numbers for each element
mod.els = [
    element_node1a(:), element_node2a(:), element_node3a(:)
    element_node1b(:), element_node2b(:), element_node3b(:)
    element_node1c(:), element_node2c(:), element_node3c(:)
    element_node1d(:), element_node2d(:), element_node3d(:)
    ];

%Now remove elements outside original boundary
[in, out] = fn_elements_in_region(mod, bdry_pts);
mod.els(out, :) = [];

%Tidy up by removing unused nodes
[mod.nds, mod.els] = fn_remove_unused_nodes(mod.nds, mod.els);

%Vectors to hold element material, element type, and absorption indices
n_els = size(mod.els, 1);
mod.el_mat_i = zeros(n_els, 1);
mod.el_typ_i = zeros(n_els, 1);
mod.el_abs_i = zeros(n_els, 1);
end