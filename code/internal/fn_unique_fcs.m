function un_fcs = fn_unique_fcs(fcs)
%Returns unique rows of fcs, taking into account that same face may be
%described by same node indices but in different order

%sort each row so that duplicate faces will be identical, but keep track of
%sort so it can be inverted after removing duplicate rows
[fcs, sort_cols] = sort(fcs, 2);

%Get the unique rows
[un_fcs, i] = unique(fcs, 'rows');
un_sort_cols = sort_cols(i, :);

%Revert the sorting on each row
[nr, nc] = size(un_fcs);
j = zeros(1, nc);
tmp = 1:nc;
for r = 1:nr
    j(un_sort_cols(r, :)) = tmp;
    un_fcs(r, :) = un_fcs(r, j);
end

end