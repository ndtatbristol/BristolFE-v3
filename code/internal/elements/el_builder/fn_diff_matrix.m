function diff_matrix = fn_diff_matrix(no_dfs, solid_or_fluid)
switch solid_or_fluid
    case 'solid'
        switch no_dfs
            case 1
                diff_matrix = [
                    1];
            case 2
                diff_matrix = [
                    1, 0
                    0, 2
                    0, 1
                    2, 0];
            case 3
                diff_matrix = [
                    1, 0, 0
                    0, 2, 0
                    0, 0, 3
                    0, 3, 2
                    3, 0, 1
                    2, 1, 0];
        end
    case 'fluid'
        diff_matrix = (1:no_dfs)';
end
end