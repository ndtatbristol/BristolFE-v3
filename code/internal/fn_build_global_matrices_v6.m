function [K, C, M, gl_lookup] = fn_build_global_matrices_v6(nds, els, el_mat_i, el_abs_i, el_typ_i, matls, el_types, fe_options)
%SUMMARY
%   Creates global matrices from mesh definitions
%INPUTS
%   nds - n x nd matrix of nodal coordinates. The row number is the node
%   number; columns are the coordinates of the node.
%   els - m x max_nds_per_el matrix of element nodes. The row number is the element
%   number; columns are the node numbers of the nodes for each element with
%   trailing zeros if element has less nodes than max_nds_per_el
%   el_mat_i - m x 1 vector of element material indices (which refer to
%   materials in 'matls' parameter
%   el_abs_i - m x 1 vector of element absorption level from 0 to 1
%   el_typ_i - m x 1 vector of element type indices
%   matls - p x 1 structured variable of materials with fields
%       matls(i).name - string giving name of material
%       matls(i).rho - density of material
%       matls(i).D - 6x6 stiffness matrix of material. 
%   [options- structured parameter of options - see default_options below for explanations]
%OUTPUTS
%   K, C, M - global (d*n) x (d*n) stiffness, damping and mass matrices, where d is no
%   of DOF at each node
%   gl_lookup - n x max(dof) matrix of global row/col indices for given
%   (node, dof) pair
%--------------------------------------------------------------------------

default_options.dof_to_use = []; %Blank uses all of available ones for all elements, subset can be set e.g. as [1,2]
%Following relate to how absorbing regions are created by adding damping
%matrix and reducing stiffness matrix to try and preserve acoustic
%impedance
default_options.damping_power_law = 3;
default_options.max_damping = 3.1415e+07; %this is the only absolute number in the list - pi * 10e6
default_options.max_stiffness_reduction = 0.01;

fe_options = fn_set_default_fields(fe_options, default_options);

%Input checks
no_els = size(els, 1);
if length(el_mat_i) ~= no_els
    error('Length of element material indices vector must equal number of elements');
end
if length(el_typ_i) ~= no_els
    error('Length of element type indices vector must equal number of elements');
end
if length(el_abs_i) ~= no_els
    error('Length of element absorbing indices vector must equal number of elements');
end

fn_console_output(sprintf('Global matrix builder v6 (nodes = %d, elements = %d, ', size(nds, 1), size(els, 1)));
t1 = clock;

%find unique element types, max DoF per element, and actual DoFs in use
unique_typ_i = unique(el_typ_i);
unique_typs = el_types(unique_typ_i);
[unique_df, max_el_df] = fn_find_dof_in_use_and_max_dof_per_el(unique_typs, fe_options.dof_to_use);
max_df = max(unique_df);

%Prepare global matrices
no_nds = size(nds, 1);
total_dof = no_nds * max_df;
K = sparse([], [], [], total_dof, total_dof);
M = sparse([], [], [], total_dof, total_dof);
C = sparse([], [], [], total_dof, total_dof);

%Find unique types, materials and shapes of elements
unique_typ_i = unique(el_typ_i); %nothing further needed as el_typ_i is the index into unique types
unique_mat_i = unique(el_mat_i); %nothing further needed as el_mat_i is the index into unique matls
[unique_el_shape_i, el_shape_i] = fn_unique_el_shape(nds, els);


%Loop over unique element types
for t = 1:numel(unique_typs)
    fn_el_mats = str2func(['fn_el_', unique_typs{t}]);
    el_i1 = el_typ_i == unique_typ_i(t); %logical indices of elements of this type

    %Find unique element matls for this type
    un_mat = unique(el_mat_i(el_i1));
    
    %Loop over unique matls for this element type
    for m = 1:numel(un_mat)
        el_i2 = (el_mat_i == un_mat(m)) & el_i1;

        if ~any(el_i2)
            %No elements of this type and material so skip to next material
            continue
        end

        %Loop over unique shapes of elements for this material and type
        for s = 1:numel(unique_el_shape_i)
            el_i3 = (el_shape_i == s) & el_i2;
            if ~any(el_i3)
                %No elements of this shape, type and material so skip to
                %next shape
                continue
            end

            %THIS IS WHERE THERE COULD BE CHUNKING LOOP?
            if un_mat(m) > 0 
                D = matls{un_mat(m)}.D;
                if isfield(matls{un_mat(m)}, 'density') %deal with legacy naming
                    rho = matls{un_mat(m)}.density;
                else
                    rho = matls{un_mat(m)}.rho;
                end
            else
                D = 0; %For elements with no material e.g. interface
                rho = 0;
            end
    
            %Get the element stiffness and mass matrices - now just for
            %single elements!
            [el_K, el_C, el_M, loc_nd, loc_df] = fn_el_mats(nds, els(unique_el_shape_i(s), :), D, rho, fe_options.dof_to_use, 0);
            %Calculate element damping matrices based on absorbing index of each element
            abs_ind = permute(el_abs_i(el_i3), [2,3,1]);
            el_M = repmat(el_M, [1, 1, numel(abs_ind)]);
            el_C = el_C + el_M .* abs_ind .^ fe_options.damping_power_law *  fe_options.max_damping;
            el_K = el_K .* exp(log(fe_options.max_stiffness_reduction) .* abs_ind .^ (fe_options.damping_power_law + 1));
            %Work out where the element matrices will go in the global matrices
            [loc_nd_i, loc_nd_j] = meshgrid(loc_nd, loc_nd);
            nd_i = reshape(els(el_i3, loc_nd_i)', size(el_K));
            nd_j = reshape(els(el_i3, loc_nd_j)', size(el_K));
            [df_i, df_j] = meshgrid(loc_df, loc_df);
            df_i = repmat(df_i, [1, 1, size(el_K, 3)]);
            df_j = repmat(df_j, [1, 1, size(el_K, 3)]);

            gi_i = (nd_i - 1) * max_df + df_i;
            gi_j = (nd_j - 1) * max_df + df_j;
            K = fn_accum_global(K, gi_i(:), gi_j(:), el_K(:));
            M = fn_accum_global(M, gi_i(:), gi_j(:), el_M(:));
            C = fn_accum_global(C, gi_i(:), gi_j(:), el_C(:));
        end
    end
end

%Reduce global matrices by knocking out any zero columns
tmp = repmat([1:no_nds], max_df, 1);
gl_nds = tmp(:);
tmp = repmat([1:max_df], 1, no_nds);
gl_dofs = tmp(:);
[gl_nds, gl_dofs, ~, K, M, C] = fn_reduce_global_matrices(gl_nds, gl_dofs, [], K, M, C);

%Produce global lookup matrix (row = node, col = DOF, content =
%global matrix index associated with node and DOF).
gl_lookup = fn_create_fast_lookup(gl_nds, gl_dofs, no_nds, 0);

ndf = size(K, 1);
tt = etime(clock, t1);
if ndf > 1e6
    fn_console_output(sprintf('DoFs = %.3fM) .......... completed in %.2f secs\n', ndf / 1e6, tt), [], 0);
else
    fn_console_output(sprintf('DoFs = %.3fk) .......... completed in %.2f secs\n', ndf / 1e3, tt), [], 0);
end

end

% function X = fn_accum_global(X, i, j, v)
% if size(v, 1) == 1
%     for j = 1:size()
% n = size(X);
% X = X + sparse(i(:), j(:), v(:), n(1), n(2));
% end

function X = fn_accum_global(X, i, j, v)
%FN_ACCUM_GLOBAL Robust sparse assembly avoiding OOM
%
%   X = fn_accum_global(X, i, j, v)
%
%   Accumulates triplet data (i,j,v) into sparse matrix X safely by
%   chunking to avoid out-of-memory errors during sparse() construction.
%
%   Inputs:
%       X : existing sparse matrix (or empty [])
%       i : row indices (vector)
%       j : column indices (vector)
%       v : values (vector)
%
%   Output:
%       X : updated sparse matrix

    % ---- parameters (tune if needed) ----
    chunkSize = 50e6;   % number of entries per chunk
    useLocalDedup = false;  % reduce duplicates per chunk

    % ---- initialise if needed ----
    if isempty(X)
        n = max(max(i), max(j));
        X = sparse(n, n);
    end

    nDOF = size(X,1);
    N = numel(v);

    for k = 1:chunkSize:N
        idx = k:min(k+chunkSize-1, N);

        ii = i(idx);
        jj = j(idx);
        vv = v(idx);

        % ---- optional: reduce duplicates within chunk ----
        if useLocalDedup
            lin = sub2ind([nDOF nDOF], ii, jj);
            [lin_u, ~, ic] = unique(lin);
            vv = accumarray(ic, vv);

            % recover i,j
            jj = ceil(lin_u / nDOF);
            ii = lin_u - (jj-1)*nDOF;
        end

        % ---- accumulate into global matrix ----
        X = X + sparse(ii, jj, vv, nDOF, nDOF);
    end
end


function [unique_el_shape_i, el_shape_i] = fn_unique_el_shape(nds,  els)
zero_vals = els == 0; %needed to avoid zero indices for elements with less than full complement of nodes
els(zero_vals) = 1;

%get flattened matrix of nodal coordinates
el_nds = reshape(nds(els(:), :), size(els, 1), size(els, 2), []);
el_nds = el_nds - el_nds(:, 1, :); %get them as relative coordinates w.r.t. first node
el_nds(repmat(zero_vals, [1, 1, size(el_nds, 3)])) = 0; %restore zeros
el_nds = reshape(el_nds, size(el_nds, 1), []);%flatten

%Find the unique rows, corresponding to the unique element shapes in the
%model and return a list of them (unique_el_shape_i) and a list of the
%associated shape index of each element in model.
[~, unique_el_shape_i, el_shape_i] = uniquetol(el_nds, 'ByRows', true);


end