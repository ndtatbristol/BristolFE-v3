function [B, detJ, N, loc_nd, loc_df] = fn_symbolic_B_and_J_matrices(shape_functions, no_df)

%should do all the shape function stuff in here and have this as the main
%entry point - inputs would be sf_powers, nodal positions in nat coords,
%no_dfs?

diff_matrix = [
    1, 0, 0
    0, 2, 0
    0, 0, 3
    0, 3, 2
    3, 0, 1
    2, 1, 0];

no_dfs = size(diff_matrix, 2);

detJ = sym('J', 'real'); %J used in file for consistency with old version but really this is detJ

no_stress = size(diff_matrix, 1);
no_J_rows = size(invJ_times_detJ, 1);
no_J_cols = size(invJ_times_detJ, 2);

el_defs = size(N, 2);

q = symvar(N);
% tmp = sym('q', [1, no_J_rows]); %name here doesn't matter - it is a dummy to ensure something to differentiate
% q = [q, tmp(no_J_rows - numel(q))];

%Calc B matrix
B = sym('B_%d%d', [no_stress, el_defs], 'real');
for i = 1:no_stress
    for j = 1:el_defs
        B(i, j) = 0;
        for k = 1:no_dfs %loop over cols in diff matrix to work out what derivatives are needed
            if diff_matrix(i, k) && diff_matrix(i, k) <= no_J_cols %means derivative w.r.t. this physical coordinate is needed
                for ii = 1:no_J_rows %loop over natural coordinate derivatives and sum after multiplying by relelvant term from inv_J
                    B(i, j) = B(i, j) + diff(N(k, j), q(ii)) * invJ_times_detJ(ii, diff_matrix(i, k)) / detJ;
                end
            end
        end
    end
end

[loc_nd, loc_df] = meshgrid([1:no_nds], [1:no_dfs]); 
loc_nd = loc_nd(:);
loc_df = loc_df(:);
end