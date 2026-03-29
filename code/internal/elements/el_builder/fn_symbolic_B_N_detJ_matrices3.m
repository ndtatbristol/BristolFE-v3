function [B, N, detJ, loc_nd, loc_df, constant_defs] = fn_symbolic_B_N_detJ_matrices3(nds_in_nat_coords, gauss_pts, gauss_weights, sf_powers, solid_or_fluid, varargin)

if numel(varargin) < 1
    simplify_expression = 0;
else
    simplify_expression = varargin{1};
end

%General constants
no_spatial_dims = 3; 
no_dims_of_el = size(gauss_pts, 2);
no_nds = size(nds_in_nat_coords, 1);
no_gps = size(gauss_pts, 1);
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
no_nds_times_dfs = no_dfs * no_nds;

%Define symbols for the physical nodal coordinates
nds = sym('nds_%d_%d', [no_nds, no_spatial_dims], 'real'); 

%Define symbols for the natural coordinates
Q = sym('q', [1, no_spatial_dims], 'real'); %nat coordinates

fprintf('Generating shape functions\n')
shape_functions = fn_symbolic_shape_functions(nds_in_nat_coords, sf_powers, Q);

N_gen = fn_symbolic_shape_function_matrix(shape_functions, no_dfs);

%Jacobian determinate and inverse Jacobian for coordinate transform
fprintf('Calculating Jacobian\n')
[detJ_gen, invJ] = fn_symbolic_inv_jacobian(shape_functions, nds, Q, no_dims_of_el);

%Calculate B-matrix (nodal displacements to strain components)
fprintf('Calculating B matrix\n')
B_gen = fn_B_matrix(N_gen, Q, solid_or_fluid, invJ);

B = sym('B_%d_%d_%d', [size(B_gen), no_gps], 'real');
N = sym('N_%d_%d_%d', [size(N_gen), no_gps], 'real');
detJ = sym('detJ_%d', no_gps, 'real');
rt3 = sym('rt3', 'real');
zero_pad = zeros(1, no_spatial_dims - no_dims_of_el);
for i = 1: no_gps
    B(:,:,i) = subs(B_gen, Q, [gauss_pts(i, :), zero_pad]);
    N(:,:,i) = subs(N_gen, Q, [gauss_pts(i, :), zero_pad]);
    detJ(i) = subs(detJ_gen, Q, [gauss_pts(i, :), zero_pad]) * gauss_weights(i);
end

%Substitute for known repeating constants
constant_defs = 'rt3 = sqrt(3);';
for i = 1: no_gps
    B(:,:,i) = subs(B(:,:,i), sqrt(3), 'rt3');
    N(:,:,i) = subs(N(:,:,i), sqrt(3), 'rt3');
    detJ(i) = subs(detJ(i), sqrt(3), 'rt3');
end

rho = sym('rho');

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

function [detJ, invJ] = fn_symbolic_inv_jacobian(shape_functions, nds, Q, no_dims_of_el)
no_dims = 3;

%Shape function matrix N for coordinates (no_dims may be less than no_dfs so local version needed here as this is only used to interpolate coordinates)
N = fn_symbolic_shape_function_matrix(shape_functions, no_dims);

%Phys coordinates in terms of natural ones
X = N * reshape(nds', [], 1);
if no_dims_of_el < no_dims
    X(no_dims_of_el + 1:end) = 0;%needs generalising for 3D
end 

J = [jacobian(X, Q), zeros(no_dims, no_dims - size(Q,2))];
if no_dims_of_el < no_dims
    J(no_dims_of_el + 1:end, no_dims_of_el + 1:end) = eye(no_dims - no_dims_of_el);%needs generalising for 3D
end 
detJ = det(J);
invJ = simplify(inv(J) * detJ) / sym('detJ', 'real');
if no_dims_of_el < no_dims
    invJ(no_dims_of_el + 1:end, :) = 0;%needs generalising for 3D
end 
end

%--------------------------------------------------------------------------

function B = fn_B_matrix(N, Q, solid_or_fluid, invJ)
diff_matrix = fn_diff_matrix(solid_or_fluid);
B = sym('B_%d%d', [size(diff_matrix, 1), size(N, 2)], 'real');
% no_dims = 3;
switch solid_or_fluid
    case 'solid'
        % no_dfs = size(diff_matrix, 2);
        for i = 1:size(diff_matrix, 1)
            for j = 1:size(N, 2)
                B(i, j) = 0;
                for k = 1:size(diff_matrix, 2) %loop over cols in diff matrix to work out what derivatives are needed
                    if diff_matrix(i, k)% && diff_matrix(i, k) <= no_Q %means derivative w.r.t. this physical coordinate is needed
                        for ii = 1:numel(Q) %loop over natural coordinate derivatives and sum after multiplying by relelvant term from inv_J
                            B(i, j) = B(i, j) + diff(N(k, j), Q(ii)) * invJ(ii, diff_matrix(i, k));
                        end
                    end
                end
            end
        end
    case 'fluid'
        N = repmat(N, [size(diff_matrix, 2), 1]);
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
    Y = Y + subs(subs(integrand, Q(1:size(gauss_pts,2)), gauss_pts(i, :)), 'detJ', detJ_i(i)) * detJ_i(i) * gauss_weights(i);
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

