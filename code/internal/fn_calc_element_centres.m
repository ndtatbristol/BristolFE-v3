function el_ctrs = fn_calc_element_centres(nds, els)
%SUMMARY
%   Returns n_els x ndims matrix of the coordinates of all element centres in a
%   model.

el_ctrs = 0;
n = 0;
for i = 1:size(els, 2) %loop of columns (node indices) of els
    j = els(:, i); %j is list of nodes referenced in this column
    v = j > 0; %v is logical vector of valid node references (i.e. non-zero and non NaN (NaN > 0 = 0) values of j)
    n = n + v; %n is vector of number of nodes per element
    j(~v) = 1; %set the invalid node indices to 1 to avoid error when referencing nodal coordinates on next line
    el_ctrs = el_ctrs + nds(j, :) .* v; %multiplication by v means that coordinates of invalid nodes are ignored anyway
end

el_ctrs = el_ctrs ./ n;

end