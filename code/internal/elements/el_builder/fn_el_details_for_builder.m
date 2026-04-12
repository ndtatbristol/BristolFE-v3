function [nds_in_nat_coords, sf_powers, gauss_pts, gauss_weights, no_dims, solid_or_fluid] = fn_el_details_for_builder(el_type, varargin)
if numel(varargin) < 1
    reduced_integration = 0;
else
    reduced_integration = varargin{1};
end
switch el_type
    case 'CPE3'
        solid_or_fluid = 'solid';
        el_shape = 'triangular';
    case 'CPE4'
        solid_or_fluid = 'solid';
        el_shape = 'quadrilateral';
    case 'C3D8'
        solid_or_fluid = 'solid';
        el_shape = 'hexahedral';
    case 'AC2D3'
        solid_or_fluid = 'fluid';
        el_shape = 'triangular';
    case 'AC2D4'
        solid_or_fluid = 'fluid';
        el_shape = 'quadrilateral';
    case 'AC3D8'
        solid_or_fluid = 'fluid';
        el_shape = 'hexahedral';
end

[nds_in_nat_coords, sf_powers, gauss_pts, gauss_weights, no_dims] = fn_el_parent_nds_and_shape_functions(el_shape, reduced_integration);

end