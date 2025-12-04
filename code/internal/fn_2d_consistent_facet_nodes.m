function [facet_nds_sorted, swapRows, feasible, singles] = fn_2d_consistent_facet_nodes(facet_nds, maxSingles)
%Following written by co-pilot with support from PW
%ORIENTPAIRS_ALLOWSINGLES  Row-wise swaps so repeated values end up in different columns,
% allowing up to 'maxSingles' values that occur only once.
%
%   [Mout, swapRows, feasible, singles] = fn_2d_consistent_facet_nodes(m)
%   [ ... ] = fn_2d_consistent_facet_nodes(m, maxSingles)
%
% INPUT:
%   m          : n-by-2 integer/numeric matrix. Each row contains a pair (a_i, b_i).
%   maxSingles : (optional) max number of single-occurrence values to allow.
%                Default = 2.
%
% OUTPUT:
%   Mout      : n-by-2 matrix after swapping some rows so that every value
%               that appears twice is in different columns across its two rows.
%               Singletons are left wherever they fall consistent with that.
%   swapRows  : logical n-by-1 indicating which rows were swapped.
%   feasible  : true if a set of swaps satisfying all "pair" constraints exists
%               with #singletons <= maxSingles; false otherwise.
%   singles   : struct describing the single-occurrence values:
%                 .values  : the singleton values (column vector)
%                 .rows    : their row indices
%                 .cols    : their original column (1 or 2)
%
% NOTES:
%   - This general formulation does NOT enforce "first occurrence in col 1".
%     It only enforces that any value with two occurrences ends up once in col 1
%     and once in col 2 in the final Mout (after swaps).
%   - Values that occur once impose no constraints; we just allow up to maxSingles of them.
%
% COMPLEXITY:
%   Near-linear in n via union-find (DSU) with parity.
%
% EXCEPTIONS / INFEASIBILITY:
%   - Any value with count > 2 ⇒ infeasible.
%   - Number of singletons > maxSingles ⇒ infeasible.
%   - Parity conflict in the constraints (rare; indicates contradictory cycle)
%     ⇒ infeasible.

    if nargin < 2 || isempty(maxSingles), maxSingles = 2; end
    if size(facet_nds,2) ~= 2
        error('Input m must be n-by-2.');
    end

    n     = size(facet_nds,1);
    vals  = facet_nds(:);                       % [m(:,1); m(:,2)] (2n x 1)
    rows  = repmat((1:n).', 2, 1);      % row index of each occurrence (2n x 1)
    colBit = [false(n,1); true(n,1)];   % col1->0, col2->1 for each occurrence

    % ---- Group by value, count occurrences
    [G, uniqueVals] = findgroups(vals);
    cnt = accumarray(G, 1);

    % Validate counts: only 1 or 2 allowed
    if any(cnt > 2)
        facet_nds_sorted = facet_nds; swapRows = false(n,1); feasible = false;
        singles = struct('values', [], 'rows', [], 'cols', []);
        return
    end

    % Identify singletons and pairs
    gSingles = find(cnt == 1);
    gPairs   = find(cnt == 2);

    if numel(gSingles) > maxSingles
        facet_nds_sorted = facet_nds; swapRows = false(n,1); feasible = false;
        % Gather singleton info to help the caller
        isSingle = ismember(G, gSingles);
        singlesIdx = find(isSingle);
        singles.values = vals(singlesIdx);
        singles.rows   = rows(singlesIdx);
        singles.cols   = 1 + colBit(singlesIdx);
        return
    end

    % Collect singleton info (handy for caller; not used in solving)
    isSingle = ismember(G, gSingles);
    singlesIdx = find(isSingle);
    singles.values = vals(singlesIdx);
    singles.rows   = rows(singlesIdx);
    singles.cols   = 1 + colBit(singlesIdx);

    % ---- Build constraints only for PAIRS (cnt == 2)
    % Use splitapply to grab the two occurrences per paired value
    occ_rows = splitapply(@(x){x}, rows,   G);
    occ_cols = splitapply(@(x){x}, colBit, G);

    % Extract (r1, c1) and (r2, c2) for the paired values
    % Each entry in gPairs has exactly two occurrences.
    r1 = zeros(numel(gPairs),1);
    r2 = zeros(numel(gPairs),1);
    c1 = false(numel(gPairs),1);
    c2 = false(numel(gPairs),1);

    for k = 1:numel(gPairs)
        rr = occ_rows{gPairs(k)};
        cc = occ_cols{gPairs(k)};
        % (rr,cc) are length 2
        r1(k) = rr(1);  r2(k) = rr(2);
        c1(k) = cc(1);  c2(k) = cc(2);
        % Optional: sanity if both occurrences land on same row (rare/odd)
        % It's harmless here (imposes no constraint) but we could check rr(1)~=rr(2)
    end

    % XOR constraint per pair: s(r1) XOR s(r2) = 1 XOR c1 XOR c2
    rhs = xor(true(size(c1)), xor(c1, c2));

    % ---- DSU with parity over rows (size n)
    parent = (1:n).';
    rank   = zeros(n,1,'uint8');
    parity = false(n,1);  % parity(i) = s(i) XOR s(parent(i))
    feasible = true;

    function [root, p2r] = dsu_find(x)
        px = parent(x);
        if px == x
            root = x; p2r = false;
        else
            [root, up] = dsu_find(px);
            parent(x) = root;
            parity(x) = xor(parity(x), up);
            p2r = parity(x);
        end
    end

    function ok = dsu_union(u, v, val)
        [ru, pu] = dsu_find(u);
        [rv, pv] = dsu_find(v);
        if ru == rv
            ok = (xor(pu, pv) == val);
            return
        end
        if rank(ru) < rank(rv)
            parent(ru) = rv;
            parity(ru) = xor(xor(pu, pv), val);
        elseif rank(ru) > rank(rv)
            parent(rv) = ru;
            parity(rv) = xor(xor(pu, pv), val);
        else
            parent(rv) = ru;
            parity(rv) = xor(xor(pu, pv), val);
            rank(ru) = rank(ru) + 1;
        end
        ok = true;
    end

    % Apply unions for all paired values; singletons impose no constraint
    for k = 1:numel(r1)
        if ~dsu_union(r1(k), r2(k), rhs(k))
            feasible = false;
            break
        end
    end

    if ~feasible
        facet_nds_sorted = facet_nds; swapRows = false(n,1);
        return
    end

    % ---- Extract swap decisions: s(i) = parity to root (take s(root)=0)
    swapRows = false(n,1);
    for i = 1:n
        [~, p] = dsu_find(i);
        swapRows(i) = p;
    end

    % ---- Apply swaps
    facet_nds_sorted = facet_nds;
    if any(swapRows)
        facet_nds_sorted(swapRows,:) = facet_nds_sorted(swapRows,[2 1]);
    end
end
