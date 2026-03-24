function [K, M, detJ, loc_nd, loc_df] = fn_symbolic_K_and_M_matrices(nds_in_nat_coords, gauss_pts, gauss_weights, sf_powers, no_dfs, solid_or_fluid)

[shape_functions, Q] = fn_symbolic_shape_functions(nds_in_nat_coords, sf_powers);

N = fn_symbolic_shape_function_matrix(shape_functions, no_dfs);

diff_matrix = fn_diff_matrix(no_dfs, solid_or_fluid);

no_dfs = size(diff_matrix, 2);
no_stress = size(diff_matrix, 1);
no_nds = size(nds_in_nat_coords, 1);
el_dfs = size(N, 2);
no_dims = size(Q, 2);

%The physical nodal coordinates
nds = sym('nds_%d_%d', [no_nds, no_dims], 'real'); 

%Jacobian determinate and inverse Jacobian for coordinate transform
[detJ_general, invJ] = fn_symbolic_inv_jacobian(shape_functions, nds, Q);

%Calculate B-matrix (nodal displacements to strain components)
B = fn_B_matrix(N, Q, diff_matrix, invJ);

%Define material stiffness matrix
D = sym('D_%d_%d', [no_stress, no_stress]); 

%Symbols for Jacobians at each Gauss point
if numel(gauss_weights) == 1
    detJ = sym('detJ_1');
else
    detJ = sym('detJ_%d', [1, numel(gauss_weights)]);
end

%Integrate to get K
K = fn_gauss_integration(B' * D * B, detJ, Q, gauss_pts, gauss_weights);

rho = sym('rho');
M = fn_gauss_integration(N' * N * rho, detJ, Q, gauss_pts, gauss_weights);
M = diag(sum(M));

%Expressions for Jacobians at each Gauss point
for i = 1:numel(gauss_weights)
    detJ(i) = simplify(subs(detJ_general, Q, gauss_pts(i, :)));
end


[loc_nd, loc_df] = meshgrid([1:no_nds], [1:no_dfs]); 
loc_nd = loc_nd(:);
loc_df = loc_df(:);
end

%--------------------------------------------------------------------------
function [shape_functions, Q] = fn_symbolic_shape_functions(nds_in_nat_coords, sf_powers)
%INPUTS
%   nds_in_nat_coords - no_nds x no_dims matrix of nodal positions in
%   natural coordinates
%   sf_powers = no_terms x no_dims matrix of powers of each natural
%   coordinate in the shape function.

no_nds = size(nds_in_nat_coords, 1);
no_dims = size(nds_in_nat_coords, 2);
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
Q = sym('q', [1, no_dims], 'real'); %nat coordinates
for n = 1:no_nds
    shape_functions(n) = 0;
    for i = 1:no_terms
        tmp = 1; %to accumlate product in the term
        for j = 1:no_dims
            tmp = tmp * Q(j) .^ sf_powers(i, j);
        end
         % * sf_coeffs(n, j)
        shape_functions(n) = shape_functions(n) + tmp * sf_coeffs(i, n);
    end
end

shape_functions = simplify(shape_functions);
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

function N = fn_symbolic_shape_function_matrix(shape_functions, no_dims)
%n is no_nodes x 1 symbolic vector of shape functions; no_dims and no_nodes
%are scalar
no_nds = numel(shape_functions);

N = sym('N', [no_dims no_dims * no_nds], 'real'); 
for i = 1:no_nds
    N(:, (i-1) * no_dims + 1: i * no_dims) = diag(repmat(shape_functions(i), [1, no_dims]));
end

end

%--------------------------------------------------------------------------

function [detJ, invJ] = fn_symbolic_inv_jacobian(shape_functions, nds, Q)
no_dims = size(nds, 2);

%Shape function matrix N for coordinates (no_dims may be less than no_dfs so local version needed here as this is only used to interpolate coordinates)
N = fn_symbolic_shape_function_matrix(shape_functions, no_dims);

%Phys coordinates in terms of natural ones
X = N * reshape(nds', [], 1);

J = jacobian(X, Q);
detJ = det(J);
invJ = inv(J);
end

%--------------------------------------------------------------------------

function B = fn_B_matrix(N, Q, diff_matrix, invJ)
no_stress = size(diff_matrix, 1);
el_dfs = size(N, 2);
no_dims = size(Q, 2);
no_dfs = size(diff_matrix, 2);
B = sym('B_%d%d', [no_stress, el_dfs], 'real');
for i = 1:no_stress
    for j = 1:el_dfs
        B(i, j) = 0;
        for k = 1:no_dfs %loop over cols in diff matrix to work out what derivatives are needed
            if diff_matrix(i, k) && diff_matrix(i, k) <= no_dims %means derivative w.r.t. this physical coordinate is needed
                for ii = 1:size(Q, 2) %loop over natural coordinate derivatives and sum after multiplying by relelvant term from inv_J
                    B(i, j) = B(i, j) + diff(N(k, j), Q(ii)) * invJ(ii, diff_matrix(i, k));
                end
            end
        end
    end
end
end

%--------------------------------------------------------------------------

function diff_matrix = fn_diff_matrix(no_dfs, solid_or_fluid)
switch solid_or_fluid
    case 'solid'
        switch no_dfs
            case 1
                diff_matrix = [
                    1];
            case 2
                diff_matrix = [
                    1, 0
                    0, 2
                    0, 1
                    2, 0];
            case 3
                diff_matrix = [
                    1, 0, 0
                    0, 2, 0
                    0, 0, 3
                    0, 3, 2
                    3, 0, 1
                    2, 1, 0];
        end
    case 'fluid'
        diff_matrix = (1:no_dfs)';
end
end

%--------------------------------------------------------------------------

function Y = fn_gauss_integration(integrand, detJ_i, Q, gauss_pts, gauss_weights)
Y = zeros(size(integrand));
for i = 1:size(gauss_pts, 1)
    Y = Y + subs(integrand, Q, gauss_pts(i, :)) * detJ_i(i) * gauss_weights(i);
end
Y = simplify(Y);
end