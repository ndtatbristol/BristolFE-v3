function [K, loc_nd, loc_df] = fn_symbolic_K_matrix(nds_in_nat_coords, gauss_pts, gauss_weights, sf_powers, no_dfs, solid_or_fluid)

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
[detJ, invJ] = fn_symbolic_inv_jacobian(shape_functions, nds, Q);

%B-matrix
B = fn_B_matrix(N, Q, diff_matrix, invJ);

%Material stiffness matrix
D = sym('D_%d_%d', [no_stress, no_stress]); 

%Integrand
integrand = simplify(B' * D * B * detJ);

%Sum over value at Gauss points
% K = sym('K_%d_%d', [el_dfs, el_dfs]); %element stiffness matrix
K = zeros(size(integrand));
for i = 1:size(gauss_pts, 1)
    K = K + subs(integrand, Q, gauss_pts(i, :)) * gauss_weights(i);
end
K = simplify(K);

[loc_nd, loc_df] = meshgrid([1:no_nds], [1:no_dfs]); 
loc_nd = loc_nd(:);
loc_df = loc_df(:);
end

function [detJ, invJ] = fn_symbolic_inv_jacobian(shape_functions, nds, Q)
no_dims = size(nds, 2);

%Shape function matrix N for coordinates (no_dims may be less than no_dfs so local version needed here as this is only used to interpolate coordinates)
N = fn_symbolic_shape_function_matrix(shape_functions, no_dims);

%Phys coordinates in terms of natural ones
X = N * reshape(nds', [], 1);

J = jacobian(X, Q);
detJ = det(J);
invJ = inv(J);
% invJ_times_detJ = invJ * detJ;
end

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