function str = fn_sprintf_numeric_matrix(X, num_fmt_str)
fmt_str = ['\t', repmat([num_fmt_str, ', '], [1, size(X, 2) - 1]), num_fmt_str, '\n'];
fmt_str = repmat(fmt_str, [1, size(X, 1)]);
str = ['[\n', sprintf(fmt_str, X), '];\n'];
end
