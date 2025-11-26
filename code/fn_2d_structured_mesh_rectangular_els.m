function mod = fn_2d_structured_mesh_rectangular_els(bdry_pts, el_size, varargin)
%USAGE
%   mod = fn_2d_structured_mesh_rectangular_els(bdry_pts, el_size)
%AUTHOR
%   Paul Wilcox (2025)
%SUMMARY
%   Utility function for generating a structured mesh of square elements
%INPUTS
%   bdry_pts - n_bdry x 2 matrix of coordinates of the n_bdry points that 
%       will define boundary of mesh.
%   el_size - element size
%   [force_square_els = 0] - force elements to be square. If bdry elements
%   define a rectangle and force_square_els = 0, element edges will be
%   exactly on boundary but elements won't be square; if force_square_els
%   = 1 then elemen edges may not be exactly on boundary but they will be
%   square.
%OUTPUT
%   mod - structured variable containing fields:
%       .nds - n_nds x 2 matrix of coordinates of each of n_nds nodes
%       .els - n_els x 4 matrix of node indices for each of n_els elements
%       .el_mat_i - n_els x 1 matrix of ones as a placeholder for element
%       material indices assigned elsewhere if more than one type of
%       material is used in model
%       .el_abs_i - n_els x 1 matrix of zeros as a placeholder for element
%       relative absorption
%       .el_typ_i - n_els x 1 matrix of zeros as a placeholder for element
%       type indices
%--------------------------------------------------------------------------
if numel(varargin) > 0
    force_square_els = varargin{1};
else
    force_square_els = 0;
end

if force_square_els
    %Figure out a bounding rectangle for whole shape that is a bit oversize
    %(excess elements will be trimmed later)
    crnr_pts = [min(bdry_pts) - 2 * el_size; max(bdry_pts) + 2 * el_size];
    x = crnr_pts(1, 1):el_size:crnr_pts(2, 1);
    y = crnr_pts(1, 2):el_size:crnr_pts(2, 2);
else
    %In this case corner points are hard on bdry
    crnr_pts = [min(bdry_pts); max(bdry_pts)];
    %Work out how many nodes are needed in x and y
    nodes_in_x_direction = ceil((crnr_pts(2, 1) - crnr_pts(1, 1)) / el_size) + 1;
    nodes_in_y_direction = ceil((crnr_pts(2, 2) - crnr_pts(1, 2)) / el_size) + 1;
    x = linspace(crnr_pts(1, 1), crnr_pts(2, 1), nodes_in_x_direction);
    y = linspace(crnr_pts(1, 2), crnr_pts(2, 2), nodes_in_y_direction);
end

%Work out nodal coordinates
[node_x_positions, node_y_positions] = meshgrid(x, y);

%Work out node numbers associated with each element
node_numbers = reshape([1:numel(node_x_positions)], numel(y), numel(x));
element_node1 = node_numbers(1:end-1, 1:end-1);
element_node2 = node_numbers(2:end,   1:end-1);
element_node3 = node_numbers(2:end,   2:end  );
element_node4 = node_numbers(1:end-1, 2:end  );

%Final m x 2 matrix of x and y coordinates for each node
mod.nds = [node_x_positions(:), node_y_positions(:)];

%Final n x 4 matrix of 4 node numbers for each element
mod.els = [element_node1(:), element_node2(:), element_node3(:), element_node4(:)];

%Now remove elements outside original boundary
[in, out] = fn_2d_find_elements_in_region(mod, bdry_pts);
mod.els(out, :) = [];

%Tidy up by removing unused nodes
[mod.nds, mod.els] = fn_remove_unused_nodes(mod.nds, mod.els);

%Associate each element with a material index = 1 and absorption index = 0 
%to start with
n_els = size(mod.els, 1);
mod.el_mat_i = ones(n_els, 1);
mod.el_abs_i = zeros(n_els, 1);
mod.el_typ_i = zeros(n_els, 1);
end