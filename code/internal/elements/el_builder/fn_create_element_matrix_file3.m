function fn_create_element_matrix_file3(fname, sym_mats)
fid = fopen(fname, 'wt');

[~, fn_name] = fileparts(fname);

%Numbers for use in this function (not directly the element script)
if isfield(sym_mats, 'B2')
    no_gps = numel(sym_mats.gauss_wts);
    K_size = size(sym_mats.B3, 2);
    D_size = size(sym_mats.D, 2);
    no_dfs = numel(unique(sym_mats.loc_df));
    no_dims = size(sym_mats.nds_sym, 2);
end
%Header line
fprintf(fid, ['function [el_K, el_C, el_M, loc_nd, loc_df] = ', fn_name, '(nds, els, D, rho, varargin)\n']);
%Comment lines
fprintf(fid, '%%SUMMARY\n');
fprintf(fid, '%%\tThis function was created automatically by fn_create_element_matrix_file3\n');
fprintf(fid, '%%\tand contains code to return the stiffness and mass matrices\n');
fprintf(fid, '%%\tfor multiple elements of the same material and type given by the latter\n');
fprintf(fid, ['%%\tpart of the filename, ', fn_name, '.\n']);
fprintf(fid, '%%INPUTS\n');
fprintf(fid, '%%\tnds - n_nds x n_dims matrix of nodal coordinates\n');
fprintf(fid, '%%\tels - n_els x n_nds_per_el matrix of node indices for each elements\n');
fprintf(fid, '%%\tD - ns x ns material stiffness matrix\n');
fprintf(fid, '%%\trho - material density\n');
fprintf(fid, '%%\t[dofs_to_use = [] - optional vector listing the DoFs to use, e.g. [1, 2]. Use [] for all]\n');
fprintf(fid, '%%OUTPUTS\n');
fprintf(fid, '%%\tel_K, el_C, el_M - n_els x n_dfs_per_el x n_dfs_per_el 3D element stiffness and mass matrices\n');
fprintf(fid, '%%AUTHOR\n');
fprintf(fid, ['%%\tPaul Wilcox (', char(datetime), ')\n']);
fprintf(fid, '\n');

%Deal with any extra starting lines specified in varargin
if isfield(sym_mats, 'start_lines')
    for i = 1:length(sym_mats.start_lines)
        fprintf(fid, [sym_mats.start_lines{i}, '\n']);
    end
end

fprintf(fid, '%%Deal with optional argument about which DOFs to use\n');
fprintf(fid, 'if isempty(varargin)\n\tdofs_to_use = [];\nelse\n\tdofs_to_use = varargin{1};\nend\n\n');

fprintf(fid, '%%Record the local node numbers of the element stiffness matrices\n');
fprintf(fid, ['loc_nd = [', num2str(sym_mats.loc_nd'), '];\n\n']);

fprintf(fid, '%%Record the local DOFs of the element stiffness matrices\n');
fprintf(fid, ['loc_df = [', num2str(sym_mats.loc_df'), '];\n\n']);

fprintf(fid, '%%Get the DOFs if not specified\n');
fprintf(fid, 'if isempty(dofs_to_use)\n\tdofs_to_use = unique(loc_df);\nend\n\n');

fprintf(fid, '%%If any inputs blank, return at this point with just the loc_nd and loc_df\n');
fprintf(fid, 'if isempty(nds) || isempty(els) || isempty(D) || isempty(rho)\n\tel_K = [];\n\tel_M = [];\n\tel_C = [];\n\t[loc_nd, loc_df] = fn_remove_dofs_from_el_matrices(loc_nd, loc_df, dofs_to_use);\n\treturn\nend\n\n');

fprintf(fid, '\n%%Some constants\n');
fprintf(fid, 'no_gauss_pts = %i;\n', no_gps);
fprintf(fid, 'no_els = size(els, 1);\n\n');

%Variables for local nds
fprintf(fid, '%%Matrices of nodal coordinates\n');
un_loc_nd = unique(sym_mats.loc_nd);
for i = 1:numel(un_loc_nd)
    for j = 1:no_dims
        fprintf(fid, [sym_mats.nds_fmt_str, ' = nds(els(:, %i), %i);\n'], un_loc_nd(i), j, un_loc_nd(i), j);
    end
end
fprintf(fid, '\n');

%Variables for Gauss weights
fprintf(fid, '%%Vector of Gauss weights\n');
fprintf(fid, 'gauss_wts = zeros(%i, 1);\n', no_gps);
for i = 1:no_gps
    fprintf(fid, 'gauss_wts(%i) = %.18e;\n', i, sym_mats.gauss_wts(i));
end

%Empty matrices for outputs
fprintf(fid, '\n%%Zero the outputs\n');
fprintf(fid, 'el_K = zeros(%i, %i, no_els);\n', K_size, K_size);
fprintf(fid, 'el_M_tmp = zeros(%i, %i, no_els);\n', K_size, K_size);
fprintf(fid, 'el_C = zeros(%i, %i, no_els);\n\n', K_size, K_size);
fprintf(fid, 'detJ = zeros(1, 1, no_els);\n');
fprintf(fid, 'B2 = zeros(%i, %i, no_els);\n', size(sym_mats.B2, 1), size(sym_mats.B2, 2));
fprintf(fid, 'B3 = zeros(%i, %i);\n', size(sym_mats.B3, 1), size(sym_mats.B3, 2));
fprintf(fid, 'N = zeros(%i, %i, no_els);\n', no_dfs, K_size);

fprintf(fid, ['%%Factors of B matrix are B1, B2, and B3. Only B2 is a function of the specific\n ' ...
    '%%element. B1 is also independent of Gauss point and is defined first.\n']);
fprintf(fid, ['B1 = ', fn_format_numeric_matrix(sym_mats.B1, '%d')]);
    
%Start loop over GPs
fprintf(fid, '%%Loop over Gauss points\n');
fprintf(fid, 'for g = 1:no_gauss_pts\n\n');

fprintf(fid, '\tswitch g\n');
for g = 1:no_gps
    fprintf(fid, '\t\tcase %i\n', g);
    fprintf(fid, '\t\t\t%%Determinant of Jacobian\n');
    fprintf(fid, fn_format_symbolic_scalar(sym_mats.detJ(g), 'detJ(:)', '\t\t\t'));

    fprintf(fid, '\t\t\t%%Terms of Jacobian\n');
    fprintf(fid, fn_format_symbolic_matrix(sym_mats.J(:, :, g), 'J', 0, '\t\t\t'));

    fprintf(fid, '\t\t\t%%Terms of B3 matrix\n');
    fprintf(fid, fn_format_symbolic_matrix(sym_mats.B3(:, :, g), 'B3', 0, '\t\t\t'));

    fprintf(fid, '\t\t\t%%Terms of N matrix\n');
    fprintf(fid, fn_format_symbolic_matrix(sym_mats.N(:, :, g), 'N', 0, '\t\t\t'));
end
fprintf(fid, '\tend\n');

fprintf(fid, '\t%%Terms of B2 matrix\n');
fprintf(fid, fn_format_symbolic_matrix(sym_mats.B2, 'B2', 0, '\t'));

fprintf(fid, '\n\t%%Calculate B matrix\n');
fprintf(fid, '\tB = pagemtimes(B1, pagemtimes(B2, B3));\n');

%Evaluate the K-matrix integrand
fprintf(fid, '\n\t%%Evaluate B''DB|J|\n');
fprintf(fid, '\tel_K = el_K + pagemtimes(pagemtimes(B, ''transpose'', pagemtimes(D, B), ''none''), detJ) * gauss_wts(g);\n\n');

%Evaluate the M-matrix integrand
fprintf(fid, '\n\t%%Evaluate rho * N''N|J|\n');
fprintf(fid, '\tel_M_tmp = el_M_tmp + rho * pagemtimes(pagemtimes(N, ''transpose'', N, ''none''), detJ) * gauss_wts(g);\n\n');

%End loop over GPs
fprintf(fid, 'end\n');

%Diagonalise M
fprintf(fid, '\n%%Diagonalise M\n');
fprintf(fid, 'el_M = zeros(%i, %i, no_els);\n', K_size, K_size);
for i = 1: K_size
    fprintf(fid, 'el_M(%i, %i, :) = sum(el_M_tmp(:, %i, :), 1);\n', i, i, i);
end

%Scale matrices (needed for fluid elements)
if isfield(sym_mats, 'scaling') && sym_mats.scaling ~= 1
    fprintf(fid, '\n%%Scale matrices\n');
    fprintf(fid, ['\nel_K = el_K * ', char(sym_mats.scaling), ';\n']);
    fprintf(fid, ['\nel_M = el_M * ', char(sym_mats.scaling), ';\n']);
end

%Call function to remove unwanted DOFs in all matrices
fprintf(fid, '\n%%Remove unwanted DOFs from element matrices\n');
fprintf(fid, 'j = ismember(loc_df, dofs_to_use);\n');
fprintf(fid, 'el_K = el_K(j, j, :);\n');
fprintf(fid, 'el_M = el_M(j, j, :);\n');
fprintf(fid, 'el_C = el_C(j, j, :);\n');
fprintf(fid, 'loc_nd = loc_nd(j);\n');
fprintf(fid, 'loc_df = loc_df(j);\n');

%Permute dimension order - this may be removed in future version but that
%requires fn_build_global_matrices to be changed
fprintf(fid, '\n%%Change dimension order of element matrices\n');
fprintf(fid, 'el_K = permute(el_K, [3, 1, 2]);\n');
fprintf(fid, 'el_M = permute(el_M, [3, 1, 2]);\n');
fprintf(fid, 'el_C = permute(el_C, [3, 1, 2]);\n');

%End function
fprintf(fid, 'end\n');

fclose(fid);
end

%--------------------------------------------------------------------------
function str = fn_format_symbolic_scalar(Z_symbolic, Z_name, varargin)
if numel(varargin) < 1
    indent_str = '';
else
    indent_str = varargin{1};
end
%This writes out symbolic matrices in flattened form ready for use in
%numeric code
if ~isscalar(Z_symbolic)
    error('Not a scalar')
end

str = [indent_str, Z_name,' = ', char(Z_symbolic), ';\n'];
str = fn_format_string_for_file(str);
end

%--------------------------------------------------------------------------
function str = fn_format_symbolic_matrix(Z_symm, var_name, symmetric, varargin)
if numel(varargin) < 1
    indent_str = '';
else
    indent_str = varargin{1};
end
%This writes out symbolic matrices in flattened form ready for use in
%numeric code
if size(Z_symm, 1) ~= size(Z_symm, 2)
    symmetric = 0;
end

str = '';

fmt_str = '(%i, %i, :)';
isvec = 0;

for i = 1:size(Z_symm,1)
    for j = 1:size(Z_symm,2)
        if ~isequal(Z_symm(i,j),sym(0))
            tmp = char(Z_symm(i,j));
            if symmetric && (j < i)
                str = [str, indent_str, sprintf([var_name, fmt_str,' = ', var_name, fmt_str, ';\n'], i, j, j, i)];
            else
                if isvec
                    str = [str, indent_str, sprintf([var_name, fmt_str,' = ', tmp, ';\n'], i)];
                else
                    str = [str, indent_str, sprintf([var_name, fmt_str,' = ', tmp, ';\n'], i, j)];
                end
            end
        end
    end
end
str = fn_format_string_for_file(str);
end

%--------------------------------------------------------------------------
function str = fn_format_string_for_file(str)
str = regexprep(str, '\^',' .^ ');
str = regexprep(str, '*',' .* ');
str = regexprep(str, '/',' ./ ');
str = regexprep(str, 'J_(\d)_(\d)', 'J($1, $2, :)');
str = [str, '\n'];
end

%--------------------------------------------------------------------------
function str = fn_format_numeric_matrix(X, num_fmt_str)
fmt_str = ['\t', repmat([num_fmt_str, ', '], [1, size(X, 2) - 1]), num_fmt_str, '\n'];
fmt_str = repmat(fmt_str, [1, size(X, 1)]);
str = ['[\n', sprintf(fmt_str, X.'), '];\n'];
end
