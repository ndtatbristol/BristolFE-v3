function L_inv = fn_invert_damping_operator(L, diag_M, tol)
% Sparse inverse of L = I + (dt/2)*inv(M)*C in O(n_s) memory via
% energy-scaled pruned Newton-Schulz, with chunked LU fallback.
% Suggested tol: 1e-12 for single precision, 1e-14 for double precision.
    ndf = size(L, 1);
    [ri, ci] = find(L);
    off = (ri ~= ci);
    if ~any(off)
        L_inv = spdiags(1 ./ full(diag(L)), 0, ndf, ndf);
        return;
    end
    coupled = unique([ri(off); ci(off)]);
    is_coupled = false(ndf, 1);
    is_coupled(coupled) = true;
    uncoupled = find(~is_coupled);
    diag_L_all = full(diag(L));
    inv_diag_unc = 1 ./ diag_L_all(uncoupled);
    L_ss   = L(coupled, coupled);
    n_s    = numel(coupled);
    sqrt_m = sqrt(abs(diag_M(coupled)));
    % 1. Primary path: energy-scaled pruned Newton-Schulz (O(n_s) memory)
    inv_L_ss = invert_schulz(L_ss, sqrt_m, tol);
    % 2. Fallback: column-chunked sparse LU (<= 64 MB dense slab)
    if isempty(inv_L_ss)
        inv_L_ss = invert_chunked_lu(L_ss, sqrt_m, tol);
    end
    [ir, ic, iv] = find(inv_L_ss);
    rows = [uncoupled; coupled(ir)];
    cols = [uncoupled; coupled(ic)];
    vals = [inv_diag_unc; iv];
    L_inv = sparse(rows, cols, vals, ndf, ndf);
end

function inv_L_ss = invert_schulz(L_ss, sqrt_m, tol)
    inv_L_ss = [];
    n_s = size(L_ss, 1);
    dL  = full(diag(L_ss));
    if any(dL == 0), return; end
    inv_dL = spdiags(1 ./ dL, 0, n_s, n_s);
    S      = spdiags(sqrt_m, 0, n_s, n_s);
    S_inv  = spdiags(1 ./ sqrt_m, 0, n_s, n_s);
    I_s    = speye(n_s);
    L_hat = S * (inv_dL * L_ss) * S_inv;
    if norm(L_hat - I_s, inf) >= 0.95
        return;
    end
    X_hat = I_s;
    two_I = 2 * I_s;
    work_tol = 0.05 * tol;
    converged = false;
    for iter = 1:16
        X_hat = prune_sparse(X_hat * (two_I - L_hat * X_hat), work_tol);
        [~, ~, rv] = find(L_hat * X_hat - I_s);
        err = max([0; abs(rv)]);
        if err <= 2 * tol
            converged = true;
            break;
        elseif err > 1
            return;
        end
    end
    if ~converged, return; end
    X_energy = prune_sparse(X_hat * inv_dL, tol);
    inv_L_ss = S_inv * X_energy * S;
end

function A = prune_sparse(A, tol)
    [r, c, v] = find(A);
    keep = (r == c) | (abs(v) >= tol);
    A = sparse(r(keep), c(keep), v(keep), size(A, 1), size(A, 2));
end


function inv_L_ss = invert_chunked_lu(L_ss, sqrt_m, tol)
    n_s = size(L_ss, 1);
    dLU = decomposition(L_ss, 'lu');  % R2017b+; or use [Lf,Uf,Pf,Qf] = lu(L_ss)
    blk = max(1, min([256, n_s, floor(64e6 / max(n_s * 8, 1))]));
    r_cell = {}; c_cell = {}; v_cell = {};
    for c0 = 1:blk:n_s
        c1 = min(c0 + blk - 1, n_s);
        b  = c1 - c0 + 1;
        rhs = zeros(n_s, b);
        rhs(sub2ind([n_s, b], c0:c1, 1:b)) = 1;
        X = dLU \ rhs;
        scaled = abs(X) .* (sqrt_m ./ sqrt_m(c0:c1)');
        mask = scaled >= tol;
        mask(sub2ind([n_s, b], c0:c1, 1:b)) = true;
        [rr, cc] = find(mask);
        r_cell{end+1} = rr; %#ok<AGROW>
        c_cell{end+1} = cc + (c0 - 1); %#ok<AGROW>
        v_cell{end+1} = X(sub2ind([n_s, b], rr, cc)); %#ok<AGROW>
    end
    inv_L_ss = sparse(vertcat(r_cell{:}), vertcat(c_cell{:}), vertcat(v_cell{:}), n_s, n_s);
end

