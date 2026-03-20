function diff_matrix = fn_diff_matrix(no_dfs, solid_or_fluid)
switch solid_or_fluid
    case 'solid'
        diff_matrix = [diag(1:no_dfs)];
    case 'fluid'
        diff_matrix = (1:no_dfs)';
end
end