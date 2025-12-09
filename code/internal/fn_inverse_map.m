function j = fn_inverse_map(i, n)
%Inverts mapping vector i so j(i) = (1:numel(i))'
%n is length of j; if not specified will be equal to max(i)
if nargin < 2
    n = max(i);
end
j = zeros(n, 1);
j(i) = 1:numel(i);
end
