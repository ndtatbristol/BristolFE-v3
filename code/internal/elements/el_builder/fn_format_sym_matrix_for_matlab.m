function str = fn_format_sym_matrix_for_matlab(Z_symm, var_name)
symmetric = 1;
%This writes out symbolic matrices in flattened form ready for use in
%numeric code
%TODO - add option to write symmetric matrix as upper triangle only with
%lower triangular elements written as copies, e.g. 
%   sprintf([var_name, '(:, %i, %i) = ', var_name, '(:, %i, %i);\n'], i, j, j, i)
str = sprintf([var_name, ' = zeros(size(els, 1), %i, %i);\n'], size(Z_symm,1), size(Z_symm,2));
for i = 1:size(Z_symm,1)
    for j = 1:size(Z_symm,2)
        if ~isequal(Z_symm(i,j),sym(0))
            tmp = char(Z_symm(i,j));
            if symmetric & (j < i)
                str = [str, sprintf([var_name, '(:, %i, %i) = ', var_name, '(:, %i, %i);\n'], i, j, j, i)];
            else
                str = [str, sprintf([var_name, '(:, %i, %i) = ', tmp, ';\n'], i, j)];
            end
        end
    end
end
str = regexprep(str, 'D_(\d)_(\d)', 'D($1, $2)');
str = regexprep(str, 'detJ_(\d)', 'detJ(:, 1, $1)');
str = regexprep(str, '\^',' .^ ');
str = regexprep(str, '*',' .* ');
str = regexprep(str, '/',' ./ ');
str = [str, '\n'];
end