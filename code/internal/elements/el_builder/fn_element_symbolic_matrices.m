function sym_mats = fn_element_symbolic_matrices(nds_in_nat_coords, gauss_pts, gauss_wts, sf_powers, solid_or_fluid, factorisation_level)
%This will be the definitive function for symbolic element matrix
%calculations. The idea is that a consistent logic is followed in all
%cases. Different factorisation levels control format of output and need to
%be implemented differently when the symbolic forms are converted into
%element files. In general, low factorisation is probably optimal for
%2D linear triangular elements, max factorisation for 3D hexahedrals.

%Factorisation levels
%   1 - return expressions for detJ, K, and M evaluated at each Gauss point
%   2 - return expressions for detJ, B, and N evaluated at each Gauss point
%   3 - return expressions for detJ, invJ, L, N, Ndiff evaluated at each Gauss point

%General constants
fprintf('Getting basics and B1 factor of B-matrix\n')
no_spatial_dims = 3; 
no_gauss_pts = size(gauss_pts, 1);
no_dims_of_el = size(gauss_pts, 2);
no_nds = size(nds_in_nat_coords, 1);
sym_mats.rho = sym('rho', 'real');
switch solid_or_fluid
    case 'solid'
        dfs = 1:3;
        no_dfs = 3;
        no_stress = 6;
        sym_mats.D = sym('D_%d_%d', [no_stress, no_stress]);
        sym_mats.B1 = [
            1 0 0  0 0 0  0 0 0
            0 0 0  0 1 0  0 0 0
            0 0 0  0 0 0  0 0 1
            0 0 0  0 0 1  0 1 0
            0 0 1  0 0 0  1 0 0
            0 1 0  1 0 0  0 0 0
            ];
        sym_mats.scaling = 1;
    case 'fluid'
        dfs = 4;
        no_stress = 1;
        sym_mats.D = sym('D'); 
        sym_mats.B1 = eye(3); %needs checking!
        sym_mats.scaling = -1 / (sym_mats.D * sym_mats.rho);
end
no_dfs = numel(dfs);
[loc_nd, loc_df] = meshgrid(1:no_nds, dfs);
sym_mats.loc_nd = loc_nd(:);
sym_mats.loc_df = loc_df(:);

%Copy Gauss weights straight to output
sym_mats.gauss_wts = gauss_wts;

%Define symbols for the physical nodal coordinates
sym_mats.nds_sym = sym('nds_%d_%d', [no_nds, no_dims_of_el], 'real'); 

%Define symbols for the natural coordinates (number = number of dimensions of element)
Q = sym('q', [1, no_dims_of_el], 'real');

fprintf('Generating shape functions\n')
shape_functions = fn_symbolic_shape_functions(nds_in_nat_coords, sf_powers, Q);

fprintf('Calculating Jacobians and determinants at Gauss points\n')
[sym_mats.detJ, sym_mats.J, sym_mats.detJ_sym, sym_mats.J_sym] = fn_jacobian_at_gauss_pts(shape_functions, sym_mats.nds_sym, Q, gauss_pts);

fprintf('Calculating B2 factor of B-matrix\n')
sym_mats.B2 = fn_calculate_B2_matrix(no_dfs, no_spatial_dims, sym_mats.J_sym, sym_mats.detJ_sym);

fprintf('Calculating B3 factor of B-matrix\n')
sym_mats.B3 = fn_calculate_B3_at_gauss_pts(shape_functions, gauss_pts, Q, no_dfs);

fprintf('Calculating N-matrices at Gauss points\n')
sym_mats.N = fn_N_at_gauss_pts(shape_functions, Q, no_dfs, gauss_pts);

if factorisation_level == 3
    return
end

sym_mats.B = sym_mats.L * kron(sym_mats.invJ, eye(3)) * kron(eye(3), sym_mats.N_diff);

%Multiply out factors to get B matrix at each Gauss point;
fprintf('Calculating B matrices at Gauss points\n')
sym_mats.B = fn_calculate_B_matrix(sym_mats.L, sym_mats.invJ, sym_mats.N_diff, no_dfs);

if factorisation_level == 2
    sym_mats = rmfield(sym_mats, {'L', 'N_diff'});
    return
end

fprintf('Calculating K matrices at Gauss points\n')
sym_mats.K = fn_calculate_K_matrix(sym_mats.B, sym_mats.D, sym_mats.detJ);
sym_mats = rmfield(sym_mats, {'B', 'N_diff'});

if nargout > 1
    test_mats = fn_form_test_mats(sym_mats, nds_in_nat_coords, eye(size(sym_mats.D)));
end
end

%--------------------------------------------------------------------------
function shape_functions = fn_symbolic_shape_functions(nds_in_nat_coords, sf_powers, Q)
%INPUTS
%   nds_in_nat_coords - no_nds x no_dims matrix of nodal positions in
%   natural coordinates
%   sf_powers = no_terms x no_dims matrix of powers of each natural
%   coordinate in the shape function.

no_nds = size(nds_in_nat_coords, 1);
no_dims_of_el = size(nds_in_nat_coords, 2);
no_dims = 3;
no_terms = size(sf_powers, 1);

%First need to get shape function coefficients, by writing simultaneous
%equations that need to be solved in order for n_i to be 1 at i^th node and
%zero at others
sf_terms = zeros(no_nds, no_nds);
for n = 1:no_nds
    sf_terms(n, :) = fn_sf_terms_at_q(sf_powers, nds_in_nat_coords(n, :));
end
sf_coeffs = inv(sf_terms);

%Create the symbolic shape functions
shape_functions = sym('n', [1, no_nds], 'real'); %shape functions
for n = 1:no_nds
    shape_functions(n) = 0;
    for i = 1:no_terms
        tmp = 1; %to accumlate product in the term
        for j = 1:no_dims_of_el
            tmp = tmp * Q(j) .^ sf_powers(i, j);
        end
        shape_functions(n) = shape_functions(n) + tmp * sf_coeffs(i, n);
    end
end

shape_functions = simplify(shape_functions);
end
%--------------------------------------------------------------------------
function [detJ, J, detJ_sym, J_sym] = fn_jacobian_at_gauss_pts(shape_functions, nds, Q, gauss_pts)

no_dims_of_el = size(Q, 2);
no_nds = numel(shape_functions);
no_gauss_pts = size(gauss_pts, 1);

%local shape function matrix just for coordinate transform to get Jacobiab
N_gen = kron(shape_functions, eye(no_dims_of_el));

%Phys coordinates in terms of natural ones
X = N_gen * reshape(nds', [], 1);

J_gen = jacobian(X, Q)';

detJ = sym('detJ_%d', [1, no_gauss_pts], 'real');
J = sym('detJ_%d', [no_dims_of_el, no_dims_of_el, no_gauss_pts], 'real');
detJ_sym = sym('detJ', 'real');
J_sym = sym('J_%d_%d', [no_dims_of_el, no_dims_of_el], 'real');

for g = 1:no_gauss_pts
    J(:, :, g) = simplify(subs(J_gen, Q, gauss_pts(g, :)));
    detJ(g) = simplify(det(J(:, :, g)));
end

end
%--------------------------------------------------------------------------
function B2 = fn_calculate_B2_matrix(no_dfs, no_spatial_dims, J_sym, detJ_sym)
no_dims_of_el = size(J_sym, 1);
B2 = kron(eye(no_dfs) ,[inv(J_sym) * det(J_sym); zeros(no_spatial_dims - no_dims_of_el, size(J_sym, 2))] / detJ_sym);
end
%--------------------------------------------------------------------------
function B3 = fn_calculate_B3_at_gauss_pts(shape_functions, gauss_pts, Q, no_dfs)
no_gauss_pts = size(gauss_pts, 1);
no_dims_of_el = size(gauss_pts, 2);
no_nds = numel(shape_functions);
%Shape function derivatives w.r.t. each dimension
% fprintf('Calculating shape function derivatives at Gauss points\n')
N_diff = fn_calculate_shape_function_derivatives_at_gauss_pts(shape_functions, gauss_pts, Q);
% fprintf('Assembling B3 factor at Gauss points\n')
B3 = sym('B3', [no_dfs * no_dims_of_el, no_dfs * no_nds, no_gauss_pts], 'real');
for g = 1:no_gauss_pts
    for n = 1: no_nds
        B3(:, (n - 1) * no_dfs + 1: n * no_dfs, g) = kron(eye(no_dfs), N_diff(:, n, g));
    end
end
end
%--------------------------------------------------------------------------
function N = fn_N_at_gauss_pts(shape_functions, Q, no_dfs, gauss_pts)
N_gen = kron(shape_functions, eye(no_dfs));
no_gauss_pts = size(gauss_pts, 1);
N = sym('N_', [size(N_gen), no_gauss_pts], 'real');
for g = 1:no_gauss_pts
    N(:, :, g) = simplify(subs(N_gen, Q, gauss_pts(g, :)));
end
end
%--------------------------------------------------------------------------
function sf_terms = fn_sf_terms_at_q(sf_powers, q)
no_dims = numel(q);
no_terms = size(sf_powers, 1);
%q(1 ... n_dim)
%sf_orders(i = 1 ... n_terms, j = 1 ... no_dims) - content indicates power
%to which component of term is raised
sf_terms = ones(1, no_terms);
for i = 1:no_terms
    for j = 1:no_dims
        sf_terms(i) = sf_terms(i) * q(j) ^ sf_powers(i, j);
    end
end
end
%--------------------------------------------------------------------------
function N_diff = fn_calculate_shape_function_derivatives_at_gauss_pts(shape_functions, gauss_pts, Q)
no_shape_functions = numel(shape_functions);
no_gauss_pts = size(gauss_pts, 1);
no_dims_of_el = size(gauss_pts, 2);
N_diff = sym('N_diff', [no_dims_of_el, no_shape_functions, no_gauss_pts], 'real');
for g = 1:no_gauss_pts
    for i = 1:no_dims_of_el
        for j = 1:no_shape_functions
            N_diff(i, j, g) = simplify(subs(diff(shape_functions(j), Q(i)), Q, gauss_pts(g, :)));
        end
    end
end
end

%--------------------------------------------------------------------------
% function E = fn_calculate_E_matrix(N_diff, no_nds, no_dfs)
% no_dims_of_el = size(N_diff, 1);
% N_diff_pad = sym('N_diff_pad', [3,3]);
% %N_diff_pad(:,:) = 0;N_diff = sym('N_diff', [no_dims_of_el, 3, 1], 'real');
% 
% N_diff_pad(1:size(N_diff, 1), :) = N_diff;
% E = kron(N_diff_pad, eye(3));
% 
% % E = sym('E', [no_dfs * no_dims_of_el, no_dfs * no_nds]);
% % E(:, :) = 0;
% % for n = 1:no_nds %outer loop across columns in steps of no_dfs
% %     for j = 1:no_dfs %outer loop down rows (first strain index - disp comp)
% %         for k = 1:no_dims_of_el %inner loop down rows (second strain index - deriv direction)
% %             E((j - 1) * no_dims_of_el + k, (n-1) * no_dfs + j) = N_diff(k, n);
% %         end
% %     end
% % end
% end
% 
%--------------------------------------------------------------------------
function B = fn_calculate_B_matrix(L, invJ, N_diff, no_dfs)


no_spatial_dims = 3;
no_gauss_pts = size(invJ, 3);
no_strain_components = size(L, 1);
no_nds = size(N_diff, 2);
% no_dims_of_el = size(N_diff, 1);
B = sym('B', [no_strain_components, no_dfs * no_nds, no_gauss_pts], 'real');
for g = 1:no_gauss_pts
    B(:, :, g) = L * kron(invJ(:,:,g), eye(no_spatial_dims)) * kron(eye(no_spatial_dims), N_diff(:,:,g));

    % % invJpadded = sym('invJpadded', [no_spatial_dims, no_dims_of_el]);
    % invJpadded = sym('invJpadded', [no_spatial_dims, no_spatial_dims]);
    % invJpadded(:, :) = 0;
    % invJpadded(1:no_dims_of_el, 1:no_dims_of_el) = invJ(:, :, g);
    % invJstar = kron(invJpadded, eye(no_dfs));
    % E = fn_calculate_E_matrix(N_diff(:, :, g), no_nds, no_dfs);
    % B(:, :, g) = L * invJstar * E;
end
end
%--------------------------------------------------------------------------

function K = fn_calculate_K_matrix(B, D, detJ)
no_gauss_pts = numel(detJ);
K = sym('B', [size(B, 2), size(B, 2), no_gauss_pts], 'real');
for g = 1:no_gauss_pts
    K(:, :, g) = simplify(B(:, :, g)' * D * B(:, :, g) * detJ(g));
end
end


