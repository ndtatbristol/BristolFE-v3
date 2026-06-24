function str = fn_format_sym_matrix_for_matlab(Z_symm, var_name, varargin)
if numel(varargin) < 1
    indent_str = '';
else
    indent_str = varargin{1};
end
symmetric = 1;
%This writes out symbolic matrices in flattened form ready for use in
%numeric code
if size(Z_symm, 1) ~= size(Z_symm, 2)
    symmetric = 0;
end

if size(Z_symm, 2) == 1
    str = sprintf([indent_str, var_name, ' = zeros(size(els, 1), %i);\n'], size(Z_symm, 1));
    fmt_str = '(:, %i)';
    symmetric = 0;
    isvec = 1;
else
    str = sprintf([indent_str, var_name, ' = zeros(size(els, 1), %i, %i);\n'], size(Z_symm,1), size(Z_symm,2));
    fmt_str = '(:, %i, %i)';
    isvec = 0;
end

for i = 1:size(Z_symm,1)
    for j = 1:size(Z_symm,2)
        if ~isequal(Z_symm(i,j),sym(0))
            tmp = char(Z_symm(i,j));
            if symmetric & (j < i)
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
str = regexprep(str, 'D_(\d)_(\d)', 'D($1, $2)');
% str = regexprep(str, 'detJ_(\d)', 'detJ(:, 1, $1)');
str = regexprep(str, 'detJ_(\d)', 'detJ(:, $1)');
str = regexprep(str, '\^',' .^ ');
str = regexprep(str, '*',' .* ');
str = regexprep(str, '/',' ./ ');
str = [str, '\n'];
end