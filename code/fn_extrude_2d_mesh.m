function mod = fn_extrude_2d_mesh(mod, z_pts)
%New function for making 3D mesh by extruding a 2D one, typically for
%waveguide problems

% Input:
%   mod.nds      : (n_nds x 2) node coordinates
%   mod.els      : (n_els x n_nds_per_el) connectivity (3 or 4)
%   mod.el_mat_i : (n_els x 1)
%   mod.el_abs_i : (n_els x 1)
%   mod.el_typ_i : (n_els x 1)
%   z_pts        : (n_z x 1) z-coordinates
%
% Output:
%   mod.nds : (n_nds*n_z x 3)
%   mod.els : (n_els*(n_z-1) x 6 or 8)

% Ensure column vector
z_pts = z_pts(:);

n_nds = size(mod.nds,1);
n_els = size(mod.els,1);
n_z   = numel(z_pts);
n_per_el = size(mod.els,2);

% ---- Create 3D nodes ----
% Repeat XY for each z level
XY = mod.nds;
mod.nds = zeros(n_nds*n_z, 3);

for k = 1:n_z
    idx = (1:n_nds) + (k-1)*n_nds;
    mod.nds(idx,1:2) = XY;
    mod.nds(idx,3)   = z_pts(k);
end

% ---- Create 3D elements ----
% Number of 3D elements
n_els_3d = n_els*(n_z-1);

if n_per_el == 3
    n_per_el_3d = 6; % wedge
elseif n_per_el == 4
    n_per_el_3d = 8; % brick
else
    error('Elements must be triangles (3) or quads (4)');
end

els3d = zeros(n_els_3d, n_per_el_3d);

% Also expand element properties
el_mat_i = repmat(mod.el_mat_i, n_z-1, 1);
el_abs_i = repmat(mod.el_abs_i, n_z-1, 1);
el_typ_i = repmat(mod.el_typ_i, n_z-1, 1);

% Build connectivity
row = 1;


for k = 1:(n_z-1)
    offset_bot = (k-1)*n_nds;
    offset_top = k*n_nds;
    
    bot = mod.els + offset_bot;
    top = mod.els + offset_top;
    
    if n_per_el == 3
        % wedge: [b1 b2 b3 t1 t2 t3]
        els3d(row:row+n_els-1, :) = [bot top];
    else
        % brick: [b1 b2 b3 b4 t1 t2 t3 t4]
        els3d(row:row+n_els-1, :) = [bot top];
    end
    
    row = row + n_els;
end

% Update mod structure
mod.els       = els3d;
mod.el_mat_i  = el_mat_i;
mod.el_abs_i  = el_abs_i;
mod.el_typ_i  = el_typ_i;

end

