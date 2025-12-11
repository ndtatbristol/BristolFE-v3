function [nds2, els2, old_nds, new_nds] = fn_sort_nds(nds1, els1, col_order)
if nargin < 3 || isempty(col_order)
    col_order = 1:size(nds1, 2);
end
[nds2, tmp] = sortrows(nds1, col_order);
new_nds = fn_inverse_map(tmp);
els2 = fn_remap_matrix(els1, new_nds);
old_nds = fn_inverse_map(new_nds);
end