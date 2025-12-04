
function T = fn_polygons_to_triangles(P, padValue)
%POLYGONS_TO_TRIANGLES Fan-triangulate rows of an m×n matrix of polygon node indices.
%
%   T = polygons_to_triangles(P)
%   T = polygons_to_triangles(P, padValue)
%
% INPUT
%   P         : m×n integer/numeric matrix. Each row lists node indices of one polygon.
%               Columns beyond the true polygon length may be padded (see padValue).
%               n must be ≥ 3; typically > 3.
%   padValue  : (optional) value used for padding. Default = 0.
%               If padValue is NaN, padding is detected via isnan().
%               If padValue is 0 (or any scalar), columns equal to padValue are ignored.
%
% OUTPUT
%   T         : K×3 integer matrix. Each row is a triangle defined by node indices.
%               K = sum_i (k_i - 2), where k_i is the number of valid vertices in row i.
%
% NOTES / ASSUMPTIONS
%   • Each row's vertices are in cyclic order around the polygon (no reordering done).
%   • We use a 'fan' from the first vertex. This partitions the polygon if the diagonals
%     from the first vertex are inside the polygon (always true for convex polygons,
%     and often acceptable for simple non-convex polygons when the first vertex is
%     a visible (convex) vertex). If you need guaranteed triangulation for arbitrary
%     non-convex polygons, you'll need coordinates and a constrained triangulation.
%
% EXAMPLE
%   % Two polygons in rows; second has padding:
%   P = [ 5  7  9 12 ;
%         3  4 10 11 14 NaN NaN ];
%   T = polygons_to_triangles(P)
%   % T will contain [5 7 9; 5 9 12; 3 4 10; 3 10 11; 3 11 14]
%
% Paul-friendly tip:
%   If your padding is zero indices, call:  T = polygons_to_triangles(P, 0);

    if nargin < 2, padValue = 0; end
    if size(P,2) < 3
        error('P must have at least 3 columns per row.');
    end

    m = size(P,1);

    % --- Detect valid vertex counts per row (exclude padding)
    if isnan(padValue)
        validMask = ~isnan(P);
    else
        validMask = P ~= padValue;
    end
    k = sum(validMask, 2);                % number of valid vertices per row

    % Rows with <3 valid vertices contribute no triangles
    kEff = max(0, k - 2);                 % triangles per row
    K = sum(kEff);                        % total triangles
    T = zeros(K, 3, class(P));            % preallocate, preserve integer class

    if K == 0
        return;
    end

    % --- Fill triangles (fan from first vertex)
    % We'll do a single pass with running write pointer.
    w = 1;
    for i = 1:m
        ki = k(i);
        if ki < 3, continue; end

        % extract this polygon's vertex list (strip padding)
        if isnan(padValue)
            verts = P(i, validMask(i,:));
        else
            verts = P(i, P(i,:) ~= padValue);
        end

        % Triangles: [v1, vj, vj+1] for j=2..ki-1
        v1 = verts(1);
        count = ki - 2;
        if count == 1
            T(w,:) = [v1, verts(2), verts(3)];
            w = w + 1;
        else
            % Vectorized fill for ki >= 4
            T(w:w+count-1, 1) = v1;
            T(w:w+count-1, 2) = verts(2:ki-1);
            T(w:w+count-1, 3) = verts(3:ki);
            w = w + count;
        end
    end
end
