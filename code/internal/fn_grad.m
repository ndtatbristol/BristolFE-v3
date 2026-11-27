function G = fn_grad(nds, els, q_nds)
%SUMMARY
%   Returns cell array of matrices G{i} (one for each dimension) that when
%   multiplied by a scalar quantity q defined at nodes q_nds gives dq/dx_i,
%   i.e. grad(q) = [G{1} * q, G{2} * q, G{3} * q] in 3D.

n_all = size(nds,1);
n_q   = numel(q_nds);
n_dim = size(nds,2);

G = cell(n_dim,1);

I = cell(n_dim,1);
J = cell(n_dim,1);
V = cell(n_dim,1);

for k = 1:n_dim
    I{k} = [];
    J{k} = [];
    V{k} = [];
end

for ii = 1:n_q

    n0 = q_nds(ii);

    nbrs = setdiff(unique(els(any(els==n0,2),:)), n0);

    A = nds(nbrs,:) - nds(n0,:);

    C = pinv(A);              % n_dim x n_nbrs

    c0 = -sum(C,2);

    for k = 1:n_dim

        % contribution from q(n0)
        I{k}(end+1) = ii;
        J{k}(end+1) = n0;
        V{k}(end+1) = c0(k);

        % contributions from neighbours
        nn = numel(nbrs);

        I{k} = [I{k} repmat(ii,1,nn)];
        J{k} = [J{k} nbrs(:).'];
        V{k} = [V{k} C(k,:)];
    end
end

for k = 1:n_dim
    G{k} = sparse(I{k},J{k},V{k},n_q,n_all);
end
end

