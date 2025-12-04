function [facet_nds, swappedRows] = fn_2d_consistent_facet_nodes(facet_nds)

%gets matrix of facet nodes into consistent ordering so that normals of
%adjacent facets are in same direction. Requires every unique node in
%facet_nds to appear exactly twice and on two different rows.

if size(facet_nds, 2) ~= 2
    error('Facet nodes must have 2 columns (line facets')
end

% function [Mout, swappedRows] = orientPairs(m)
n    = size(facet_nds,1);
vals = facet_nds(:);                             % column-major: [m(:,1); m(:,2)]
rows = repelem((1:n).', 2);              % row index for each occurrence
% cols = repmat([1; 2], n, 1);          % (not needed explicitly)

% ---- Group by value, check "appears exactly twice"
[G, ~] = findgroups(vals);
cnt = accumarray(G, 1);
if ~all(cnt == 2)
    error('Each unique value must appear exactly twice on two different rows.');
end

% ---- For each value, find its first (earlier) row
firstRowPerVal = splitapply(@min, rows, G);   % one per group/value

% Desired column for each occurrence: 1 if it is at the first row, else 2
desiredCol = 1 + (rows ~= firstRowPerVal(G));  % vector of 1s/2s, length 2n

% Reshape to n-by-2 so columns align with current m(:,1) and m(:,2)
desiredCol = reshape(desiredCol, n, 2);

% Rows already correct ( [1 2] ), rows needing swap ( [2 1] )
okRows    =  (desiredCol(:,1) == 1) & (desiredCol(:,2) == 2);
swapRows  =  (desiredCol(:,1) == 2) & (desiredCol(:,2) == 1);

% Any impossible rows? (both firsts => [1 1], both seconds => [2 2])
badRows = ~(okRows | swapRows);
if any(badRows)
    r = find(badRows, 1, 'first');
    error('Row %d is inconsistent (both entries first or both second). Cannot satisfy requirement.', r);
end

% ---- Apply swaps
facet_nds(swapRows, :) = facet_nds(swapRows, [2 1]);

if nargout > 1
    swappedRows = swapRows;
end
end


