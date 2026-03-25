function N = fn_symbolic_shape_function_matrix(shape_functions, no_dims)
%n is no_nodes x 1 symbolic vector of shape functions; no_dims and no_nodes
%are scalar
no_nds = numel(shape_functions);

N = sym('N', [no_dims no_dims * no_nds], 'real'); 
for i = 1:no_nds
    N(:, (i-1) * no_dims + 1: i * no_dims) = diag(repmat(shape_functions(i), [1, no_dims]));
end

end