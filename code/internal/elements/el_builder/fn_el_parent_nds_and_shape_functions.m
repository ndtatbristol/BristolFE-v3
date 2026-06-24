function [nds_in_nat_coords, sf_powers, gauss_pts, gauss_weights, no_dims] = fn_el_parent_nds_and_shape_functions(shape, varargin)
if numel(varargin) < 1
    reduced_integration = 0;
else
    reduced_integration = varargin{1};
end

switch shape
    %1D elements
    case 'line'
        nds_in_nat_coords = [
            0
            1];
        sf_powers = [
            0
            1];
        gauss_pts = 1/2;
        gauss_weights = 1;
    %2D elements
    case 'triangular'
        nds_in_nat_coords = [
            0, 0
            1, 0
            0, 1];
        sf_powers = [
            0, 0
            1, 0
            0, 1];
        gauss_pts = [
            1, 1] / 3;
        gauss_weights = 1/2;
    case 'quadrilateral'
        nds_in_nat_coords = [
            -1, -1
            -1,  1
            1,  1
            1, -1];
        sf_powers = [
            0, 0
            1, 0
            0, 1
            1, 1];
        if reduced_integration
            gauss_pts = [
                0, 0];
            gauss_weights = 4;
        else
            gauss_pts = nds_in_nat_coords / sqrt(3);
            gauss_weights = ones(size(gauss_pts, 1), 1);
        end
%3D elements
    case 'tetrahedral'
        nds_in_nat_coords = [
            0, 0, 0
            1, 0, 0
            0, 1, 0
            0, 0, 1];
        sf_powers = [
            0, 0, 0
            1, 0, 0
            0, 1, 0
            0, 0, 1];
        gauss_pts = [
            1, 1, 1] / 4;
        gauss_weights = 1/6;

    case 'hexahedral'
        nds_in_nat_coords = [
            -1  -1  -1
            1  -1  -1
            1   1  -1
            -1   1  -1
            -1  -1   1
            1  -1   1
            1   1   1
            -1   1   1
            ];
        sf_powers = [
            0, 0, 0
            1, 0, 0
            0, 1, 0
            0, 0, 1
            0, 1, 1
            1, 0, 1
            1, 1, 0
            1, 1, 1];
        if reduced_integration
            gauss_pts = [
                0, 0, 0];
            gauss_weights = 8;
        else
            gauss_pts = nds_in_nat_coords / sqrt(3);
            gauss_weights = ones(size(gauss_pts, 1), 1);
        end
    case 'triangular_prism'
        nds_in_nat_coords = [
            0   0  -1
            1   0  -1
            0   1  -1
            0   0   1
            1   0   1
            0   1   1
            ];
        sf_powers = [
            0, 0, 0
            1, 0, 0
            0, 1, 0
            0, 0, 1
            0, 1, 1
            1, 0, 1
            ];
        gauss_pts = [
            1/6  1/6  -1/sqrt(3)
            2/3  1/6  -1/sqrt(3)
            1/6  2/3  -1/sqrt(3)
            1/6  1/6   1/sqrt(3)
            2/3  1/6   1/sqrt(3)
            1/6  2/3   1/sqrt(3)
            ];
        gauss_weights = ones(size(gauss_pts, 1), 1) / 6;
end
no_dims = size(nds_in_nat_coords, 2);

end