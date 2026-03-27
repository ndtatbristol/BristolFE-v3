function [B, D, detJ_general, Q, N, loc_nd, loc_df, constant_defs] = fn_symbolic_B_D_J_Q_matrices(nds_in_nat_coords, gauss_pts, gauss_weights, sf_powers, solid_or_fluid, varargin)

if numel(varargin) < 1
    simplify_expression = 0;
else
    simplify_expression = varargin{1};
end

%Define material stiffness matrix
switch solid_or_fluid
    case 'solid'
        no_dfs = 3;
        no_stress = 6;
        D = sym('D_%d_%d', [no_stress, no_stress]); 
    case 'fluid'
        no_dfs = 1;
        fluid_dof = 4;
        no_stress = 1;
        D = sym('D'); 
end

fprintf('Generating shape functions\n')
[shape_functions, Q] = fn_symbolic_shape_functions(nds_in_nat_coords, sf_powers);

% candidates = fn_factor_shape_function(shape_functions(1));

no_spatial_dims = 3; %size(Q, 2);
no_dims_of_el = size(Q, 2);
N = fn_symbolic_shape_function_matrix(shape_functions, no_spatial_dims);
diff_matrix = fn_diff_matrix(solid_or_fluid);

no_dfs = size(diff_matrix, 2); %OK
no_nds = size(nds_in_nat_coords, 1); %OK
no_nds_times_dfs = no_dfs * no_nds;

%The physical nodal coordinates
nds = sym('nds_%d_%d', [no_nds, no_dims_of_el], 'real'); 

%Jacobian determinate and inverse Jacobian for coordinate transform
fprintf('Calculating Jacobian\n')
[detJ_general, invJ] = fn_symbolic_inv_jacobian(shape_functions, nds, Q);

%Calculate B-matrix (nodal displacements to strain components)
fprintf('Calculating B matrix\n')
B = fn_B_matrix(N, Q, diff_matrix, invJ, no_nds_times_dfs);

constant_defs = 'root3 = sqrt(3);';
root3 = sym('root3', 'real');
for i = 1:numel(B)
    B(i) = subs(B(i), sqrt(3), 'root3');
end

%Symbols for Jacobians at each Gauss point
if isscalar(gauss_weights)
    detJ = sym('detJ_1');
else
    detJ = sym('detJ_%d', [1, numel(gauss_weights)]);
end


fprintf('Calculating M matrix\n')
rho = sym('rho');
M = fn_gauss_integration(N' * N * rho, detJ, Q, gauss_pts, gauss_weights, simplify_expression);
M = diag(sum(M));

% if strcmp(solid_or_fluid, 'fluid')
%     K = -K / rho / D;
%     M = -M / D / rho;
% end
% 
% %Expressions for Jacobians at each Gauss point
% fprintf('Calculating Jacobians at Gauss points\n')
% detJ = fn_jacobians_at_gauss_points(detJ_general, Q, gauss_pts);

switch solid_or_fluid
    case 'solid'
        [loc_nd, loc_df] = meshgrid(1:no_nds, 1:no_dfs); 
    case 'fluid'
        [loc_nd, loc_df] = meshgrid([1:no_nds], fluid_dof); 
end
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
invJ = simplify(inv(J) * detJ) / sym('detJ', 'real');
 
end

%--------------------------------------------------------------------------

function B = fn_B_matrix(N, Q, diff_matrix, invJ, no_nds_times_dfs)
no_stress = size(diff_matrix, 1);
no_dims = size(Q, 2);
no_dfs = size(diff_matrix, 2);
B = sym('B_%d%d', [no_stress, no_nds_times_dfs], 'real');
for i = 1:no_stress
    for j = 1:no_nds_times_dfs
        B(i, j) = 0;
        for k = 1:no_dims %loop over cols in diff matrix to work out what derivatives are needed
            if diff_matrix(i, k) && diff_matrix(i, k) <= no_dims %means derivative w.r.t. this physical coordinate is needed
                for ii = 1:no_dims %loop over natural coordinate derivatives and sum after multiplying by relelvant term from inv_J
                    B(i, j) = B(i, j) + diff(N(k, j), Q(ii)) * invJ(ii, diff_matrix(i, k));
                end
            end
        end
    end
end
end

%--------------------------------------------------------------------------

function diff_matrix = fn_diff_matrix(solid_or_fluid)
switch solid_or_fluid
    case 'solid'
        diff_matrix = [
            1, 0, 0
            0, 2, 0
            0, 0, 3
            0, 3, 2
            3, 0, 1
            2, 1, 0];
    case 'fluid'
        diff_matrix = [1, 2, 3];
end
end

%--------------------------------------------------------------------------

function Y = fn_gauss_integration(integrand, detJ_i, Q, gauss_pts, gauss_weights, varargin)
if ~isempty(varargin)
    simplify_expr = varargin{1};
else
    simplify_expr = false;
end
Y = zeros(size(integrand));
for i = 1:size(gauss_pts, 1)
    Y = Y + subs(subs(integrand, Q, gauss_pts(i, :)), 'detJ', detJ_i(i)) * detJ_i(i) * gauss_weights(i);
    fprintf('  integrating %i of %i\n', i, size(gauss_pts, 1));
end
if simplify_expr
    for i = 1:numel(Y)
        Y(i) = simplify(Y(i), simplify_expr);
        fprintf('  simplifying %i of %i\n', i, numel(Y));
    end
end
end

function Y = fn_gauss_integration2(integrand_terms, detJ_i, Q, gauss_pts, gauss_weights, varargin)
if ~isempty(varargin)
    simplify_expr = varargin{1};
else
    simplify_expr = false;
end
Y = zeros(size(integrand_terms{1}, 1), size(integrand_terms{end}, 2));
for i = 1:size(gauss_pts, 1)
    tmp = 1;
    for j = 1:numel(integrand_terms)
        tmp = tmp * simplify(subs(integrand_terms{j}, Q, gauss_pts(i, :)));
    end
    Y = Y + subs(tmp, 'detJ', detJ_i(i)) * detJ_i(i) * gauss_weights(i);
    fprintf('  integrating %i of %i\n', i, size(gauss_pts, 1));
end
if simplify_expr
    for i = 1:numel(Y)
        Y(i) = simplify(Y(i), simplify_expr);
        fprintf('  simplifying %i of %i\n', i, numel(Y));
    end
end
end

%--------------------------------------------------------------------------
function detJ = fn_jacobians_at_gauss_points(detJ_general, Q, gauss_pts)
detJ = sym('detJ', [size(gauss_pts, 1), 1], 'real')
for i = 1:size(gauss_pts, 1)
    detJ(i) = simplify(subs(detJ_general, Q, gauss_pts(i, :)));
end
end

