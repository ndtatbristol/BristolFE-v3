function fn_create_element_matrix_file2(fname, B, D, J, N, Q, gauss_pts, gauss_weights, loc_nd, loc_df, no_dims, varargin)
fid = fopen(fname, 'wt');

[~, fn_name] = fileparts(fname);

%Header line
fprintf(fid, ['function [el_K, el_C, el_M, loc_nd, loc_df] = ', fn_name, '(nds, els, D, rho, varargin)\n']);
%Comment lines
fprintf(fid, '%%SUMMARY\n');
fprintf(fid, '%%\tThis function was created automatically by fn_create_element_matrix_file\n');
fprintf(fid, '%%\tand contains code to return the stiffness and mass matrices\n');
fprintf(fid, '%%\tfor multiple elements of the same material and type given by the latter\n');
fprintf(fid, ['%%\tpart of the filename, ', fn_name, '.\n']);
fprintf(fid, '%%INPUTS\n');
fprintf(fid, '%%\tnds - n_nds x n_dims matrix of nodal coordinates\n');
fprintf(fid, '%%\tels - n_els x n_nds_per_el matrix of node indices for each elements\n');
fprintf(fid, '%%\tD - ns x ns material stiffness matrix\n');
fprintf(fid, '%%\trho - material density\n');
fprintf(fid, '%%\t[dofs_to_use = [] - optional string listing the DoFs to use, e.g. ''12''. Use [] for all]\n');
fprintf(fid, '%%OUTPUTS\n');
fprintf(fid, '%%\tel_K, el_C, el_M - n_els x n_dfs_per_el x n_dfs_per_el 3D element stiffness and mass matrices\n');
fprintf(fid, '%%AUTHOR\n');
fprintf(fid, ['%%\tPaul Wilcox (', char(datetime), ')\n']);
fprintf(fid, '\n');

%Deal with any extra lines specified in varargin
for i = 1:length(varargin)
    fprintf(fid, [varargin{i}, '\n']);
end

fprintf(fid, '%%Deal with optional argument about which DOFs to use\n');
fprintf(fid, 'if isempty(varargin)\n\tdofs_to_use = [];\nelse\n\tdofs_to_use = varargin{1};\nend\n\n');

fprintf(fid, '%%Record the local node numbers of the element stiffness matrices\n');
fprintf(fid, ['loc_nd = [', num2str(loc_nd'), '];\n\n']);

fprintf(fid, '%%Record the local DOFs of the element stiffness matrices\n');
fprintf(fid, ['loc_df = [', num2str(loc_df'), '];\n\n']);

fprintf(fid, '%%Get the DOFs if not specified\n');
fprintf(fid, 'if isempty(dofs_to_use)\n\tdofs_to_use = unique(loc_df);\nend\n\n');

fprintf(fid, '%%If any inputs blank, return at this point with just the loc_nd and loc_df\n');
fprintf(fid, 'if isempty(nds) || isempty(els) || isempty(D) || isempty(rho)\n\tel_K = [];\n\tel_M = [];\n\tel_C = [];\n\tel_Q = [];\n\t[loc_nd, loc_df] = fn_remove_dofs_from_el_matrices(loc_nd, loc_df, dofs_to_use);\n\treturn\nend\n\n');

%Shortcut variables for local nds
fprintf(fid, '%%Temporary matrices of nodal coordinates to save time\n');
un_loc_nd = unique(loc_nd);
for i = 1:numel(un_loc_nd)
    for j = 1:no_dims
        fprintf(fid, 'nds_%i_%i = nds(els(:, %i), %i);\n', un_loc_nd(i), j, un_loc_nd(i), j);
    end
end
fprintf(fid, '\n');

%Start loop over Gauss points
K_size = size(B, 2);
D_size = size(D, 1);
fmt_str = [repmat('%.10g ', 1, size(gauss_pts, 2))];
fprintf(fid, '\n%%Define Gauss points and weights\n');
fprintf(fid,     'Q = [');
for i = 1:size(gauss_pts, 1)
    fprintf(fid, fmt_str, gauss_pts(i, :));
    if i < size(gauss_pts, 1)
        fprintf(fid, '\n    ');
    else
        fprintf(fid, '];\n');
    end
end
fprintf(fid,     'W = [');
for i = 1:size(gauss_pts, 1)
    fprintf(fid, '%.10g', gauss_weights(i));
    if i < size(gauss_pts, 1)
        fprintf(fid, '\n    ');
    else
        fprintf(fid, '];\n');
    end
end

fprintf(fid, '\n%%Zero the outputs\n');
fprintf(fid, 'el_K = zeros(size(els, 1), %i, %i);\n', K_size, K_size);
fprintf(fid, 'tmp_M = zeros(size(els, 1), %i, %i);\n', K_size, K_size);
fprintf(fid, 'el_C = zeros(size(els, 1), %i, %i);\n', K_size, K_size);
fprintf(fid, 'BTDB = zeros(size(els, 1), %i, %i);\n', K_size, K_size);
fprintf(fid, 'BTD = zeros(size(els, 1), %i, %i);\n', K_size, D_size);
fprintf(fid, 'NTN = zeros(size(els, 1), %i, %i);\n', K_size, K_size);

fprintf(fid, '\n%%Loop over Gauss points\n');
fprintf(fid, 'for i = 1:size(Q, 1)\n');
%Define the Gauss point variables q1, q2, ... used in equations
for j = 1:size(gauss_pts, 2)
    fprintf(fid, '    q%i = Q(i, %i);\n', j, j);
end

%Evaluate Jacobian
fprintf(fid, '\n    %%Evaluate Jacobian\n');
fprintf(fid, ['    detJ = ', regexprep(char(J), '*',' .* '), ';\n']);

%Evaluate terms in B matrix
fprintf(fid, '\n    %%Evaluate B matrix\n');
fprintf(fid, fn_format_sym_matrix_for_matlab(B, 'B', '    '));

%Get symbolic expression for B'*D
fprintf(fid, '\n    %%Evaluate B''D\n');
for i = 1:K_size
    for j = 1: D_size
        fprintf(fid, '    BTD(:, %i, %i) = ', i, j);
        for k = 1: D_size
            fprintf(fid, ' + B(:, %i, %i) * D(%i, %i)', k, i, k, j);
        end
        fprintf(fid, ';\n');
    end
end

%Get symbolic expression for B'*D*B
fprintf(fid, '\n    %%Evaluate B''DB\n');
for i = 1:K_size
    for j = i: K_size
        fprintf(fid, '    BTDB(:, %i, %i) = ', i, j);
        for k = 1: D_size
            fprintf(fid, ' + BTD(:, %i, %i) .* B(:, %i, %i)', i, k, k, j);
        end
        fprintf(fid, ';\n');
    end
end

%Multiply by Jacobian and Gauss weight and add to K
fprintf(fid, '\n    %%Evaluate contribution to K at Gauss point and accumulate\n');
for i = 1:K_size
    for j = i: K_size
        fprintf(fid, '    el_K(:, %i, %i) = el_K(:, %i, %i) + BTDB(:, %i, %i) .* detJ * W(i);\n', i, j, i, j, i, j);
    end
end

%Evaluate terms in N matrix
fprintf(fid, '\n    %%Evaluate N matrix\n');
fprintf(fid, fn_format_sym_matrix_for_matlab(N, 'N', '    '));

%Get symbolic expression for N'*N
fprintf(fid, '\n    %%Evaluate N''N\n');
for i = 1: K_size
    for j = i: K_size
        fprintf(fid, '    NTN(:, %i, %i) = ', i, j);
        for k = 1: size(N, 1)
            fprintf(fid, ' + N(:, %i, %i) .* N(:, %i, %i)', k, i, k, j);
        end
        fprintf(fid, ';\n');
    end
end

fprintf(fid, '\n    %%Evaluate contribution to M at Gauss point and accumulate\n');
for i = 1:K_size
    for j = i: K_size
        fprintf(fid, '    tmp_M(:, %i, %i) = tmp_M(:, %i, %i) + NTN(:, %i, %i) .* rho * W(i);\n', i, j, i, j, i, j);
    end
end

%End of Gauss point loop
fprintf(fid, 'end\n\n');

fprintf(fid, '\n%%Copy upper triangles of K and tmp_M into lower for symmmetry\n');
for i = 1: K_size
    for j = 1: i - 1
        fprintf(fid, 'el_K(:, %i, %i) = el_K(:, %i, %i);\n', i, j, j, i);
        fprintf(fid, 'tmp_M(:, %i, %i) = tmp_M(:, %i, %i);\n', i, j, j, i);
    end
end

fprintf(fid, '\n%%Diagonalise M\n');
fprintf(fid, 'M = zeros(size(els, 1), %i, %i);\n', K_size, K_size);
for i = 1: K_size
    fprintf(fid, 'el_M(:, %i, %i) = sum(tmp_M(:, :, %i));\n', i, i, i);
    % for j = 1: K_size
    %     fprintf(fid, ' + tmp_M(:, %i, %i)', j, i);
    % end
    % fprintf(fid, ';\n');
end

%Call function to remove unwanted DOFs in all matrices
fprintf(fid, '%%CRemove unwanted DOFs from element matrices\n');
fprintf(fid, '[loc_nd, loc_df, el_K, el_C, el_M] = fn_remove_dofs_from_el_matrices(loc_nd, loc_df, dofs_to_use, el_K, el_C, el_M);\n');

%End line and close
fprintf(fid, '\nend\n');
fclose(fid);


end