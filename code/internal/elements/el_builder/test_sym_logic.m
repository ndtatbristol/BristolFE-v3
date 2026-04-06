clear
no_nds = 4;
no_dfs = 3;
no_dims_of_el = 2;
no_spatial_dims = 3;
Nd = sym('dN%d_by_dq%d', [no_nds, no_dims_of_el],  'real')';
v = sym('v', [no_nds, no_dfs],  'real')';
v = v(:);
B3 = sym('B3', [no_dfs * no_dims_of_el, no_dfs * no_nds]);
for n = 1:no_nds
    B3(:, (n - 1) * no_dfs + 1: n * no_dfs) = kron(eye(no_dfs), Nd(:, n));
end
% B3 = kron( Nd, eye(no_dfs));%this is wrong!!!
invJ = sym('dq%d_by_dx%d', [no_dims_of_el, no_dims_of_el],  'real')';
invJpad = sym('invJ', [no_spatial_dims, no_dims_of_el],  'real');
invJpad(:, :) = 0;
invJpad(1:no_dims_of_el, :) = invJ;
B2 = kron(eye(no_spatial_dims), invJpad);
B1= [
    1 0 0  0 0 0  0 0 0
    0 0 0  0 1 0  0 0 0
    0 0 0  0 0 0  0 0 1
    0 0 0  0 0 1  0 1 0
    0 0 1  0 0 0  1 0 0
    0 1 0  1 0 0  0 0 0
    ];

 B = B1 * B2 * B3;%this should be it!