function mod = fn_2d_structured_mesh(bdry_pts, el_size, el_type, varargin)
%USAGE
%   mod = fn_2d_structured_mesh(bdry_pts, el_size, el_type [, force_equilateral_els])
%AUTHOR
%   Paul Wilcox (2026)
%SUMMARY
%   Utility function for generating a structured mesh of triangular or
%   quadrilateral elements, that fills the region specified by bdry_nds.
%INPUTS
%   bdry_pts - n_bdry x 2 matrix of coordinates of the n_bdry points that 
%       will define boundary of mesh.
%   el_size - element size
%   el_type - valid 2D element type. Number of nodes of element type will
%   detemine if mesh is for triangualar or quadrilateral elements.
%   [force_equilateral_els = 0] - force elements to be equilaterial 
%   triangles or squares. If bdry elements define a rectangle and 
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
tmp = fn_query_el_type_info(el_type);
switch tmp.shape
    case 'triangular'
        mod = fn_2d_structured_mesh_triangular_els(bdry_pts, el_size, el_type, varargin{:});
    case 'quadrilateral'
        mod = fn_2d_structured_mesh_rectangular_els(bdry_pts, el_size, el_type, varargin{:});
    otherwise
        error([el_type, ' is unknown or not a 2D element type']);
end

end