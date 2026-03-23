function [shape_functions, Q] = fn_shape_functions_sym(nds_in_nat_coords, sf_powers)
%INPUTS
%   nds_in_nat_coords - no_nds x no_dims matrix of nodal positions in
%   natural coordinates
%   sf_powers = no_terms x no_dims matrix of powers of each natural
%   coordinate in the shape function.

no_nds = size(nds_in_nat_coords, 1);
no_dims = size(nds_in_nat_coords, 2);
no_terms = size(sf_powers, 1);
% no_dims = 3;
% if size(sf_powers, 2) < no_dims
%     sf_powers = [sf_powers, zeros(size(sf_powers, 1), no_dims - size(sf_powers, 2))];
% end

%First need to get shape function coefficients, by writing simultaneous
%equations that need to be solved in order for n_i to be 1 at i^th node and
%zero at others
sf_terms = zeros(no_nds, no_nds);
for n = 1:no_nds
    sf_terms(n, :) = fn_sf_terms_at_q(sf_powers, nds_in_nat_coords(n, :));
end
sf_coeffs = inv(sf_terms);

%Create the symbolic shape functions
shape_functions = sym('n', [1, no_nds], 'real'); %shape functions
Q = sym('q', [1, no_dims], 'real'); %nat coordinates
for n = 1:no_nds
    shape_functions(n) = 0;
    for i = 1:no_terms
        tmp = 1; %to accumlate product in the term
        for j = 1:no_dims
            tmp = tmp * Q(j) .^ sf_powers(i, j);
        end
         % * sf_coeffs(n, j)
        shape_functions(n) = shape_functions(n) + tmp * sf_coeffs(i, n);
    end
end

shape_functions = simplify(shape_functions);
end

function sf_terms = fn_sf_terms_at_q(sf_powers, q)
no_dims = numel(q);
no_terms = size(sf_powers, 1);
%q(1 ... n_dim)
%sf_orders(i = 1 ... n_terms, j = 1 ... no_dims) - content indicates power
%to which component of term is raised
sf_terms = ones(1, no_terms);
for i = 1:no_terms
    for j = 1:no_dims
        sf_terms(i) = sf_terms(i) * q(j) ^ sf_powers(i, j);
    end
end
end
