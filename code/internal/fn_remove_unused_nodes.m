function [nds2, els2, old_nds, new_nds] = fn_remove_unused_nodes(nds, els)

%Note how the x-refs work
%   nds2 = nds(old_nds, :);
%   els2 = reshape(new_nds(els), size(els));
%Conceptually:
%   nds = nds2(new_nds, :); %except this actual expression won't work 
%   because the entries corresponding to deleted nodes in
%   new_nodes have zero indices. However, it can be used to remap anything
%   else that references only the old node numbers:
%   X2 = reshape(new_nds(X), size(X));

%list of original node numbers
new_nds = [1:size(nds, 1)]';
old_nds = [1:size(nds, 1)]';
in_use = zeros(size(old_nds));

u = unique(els(:));
u(u == 0 | isnan(u)) = [];
in_use(u) = 1;
in_use = find(in_use);
nds2 = nds(in_use, :);
old_nds = old_nds(in_use);

new_nds = fn_inverse_map(old_nds, size(nds, 1));
els2 = fn_remap_matrix(els, new_nds);

end