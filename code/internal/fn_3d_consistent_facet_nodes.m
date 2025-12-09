function [facet_nds_out, eds, ed_fcs] = fn_3d_consistent_facet_nodes(facet_nds)
%Written by co-pilot

% function facet_nds_out = fn_3d_consistent_facet_nodes(facet_nds)
%FN_3D_CONSISTENT_FACET_NODES Reorder triangle nodes to enforce opposite orientation
% across shared edges.
%
%   facet_nds_out = fn_3d_consistent_facet_nodes(facet_nds)
%
% INPUT:
%   facet_nds : (n x 3) integer array. Each row [A B C] lists the node indices
%               for one triangular facet. The three consecutive pairs (A,B),
%               (B,C), (C,A) are the oriented edges of the triangle. Each
%               undirected edge {i,j} appears in at most two rows.
%
% OUTPUT:
%   facet_nds_out : (n x 3) integer array (same triangles) but with rows
%                   possibly reordered so that any shared edge between two
%                   triangles appears with reversed node order (i.e., if one
%                   triangle lists edge u->v consecutively, the other lists v->u).
%
% NOTES:
% - Boundary edges (edges appearing only once) are left as-is.
% - If the mesh is non-orientable or contains a conflict (e.g., cycles that
%   force contradictory orientations), the function will issue a warning for the
%   conflicting edge/triangle but will otherwise leave previously oriented
%   triangles unchanged.
%
% Example:
%   facet_nds = [1 2 3; 2 3 4];
%   % Triangles share edge {2,3}; both initially have 2->3.
%   out = fn_3d_consistent_facet_nodes(facet_nds)
%   % out = [1 2 3; 3 2 4]  (now the shared edge is reversed in the 2nd row)
%

    % -------------------- Input checks --------------------
    if ~isnumeric(facet_nds) || size(facet_nds,2) ~= 3
        error('facet_nds must be an (n x 3) numeric matrix.');
    end
    if any(~isfinite(facet_nds(:))) || any(facet_nds(:) ~= round(facet_nds(:)))
        error('facet_nds must contain integer node indices.');
    end

    n = size(facet_nds,1);
    facet_nds_out = facet_nds;

    % -------------------- Build edge->triangle map --------------------
    % Map key is 'min-max' (string), value is array of triangle indices {i,j}
    edgeMap = containers.Map('KeyType','char','ValueType','any');

    for i = 1:n
        tri = facet_nds_out(i,:);
        edges = [tri([1 2]); tri([2 3]); tri([3 1])]; % 3x2: oriented edges
        for k = 1:3
            u = edges(k,1); v = edges(k,2);
            key = edge_key(u, v);
            if isKey(edgeMap, key)
                edgeMap(key) = [edgeMap(key), i];
            else
                edgeMap(key) = i;
            end
        end
    end

    % -------------------- BFS to propagate consistent orientation --------------------
    visited = false(n,1);

    for start = 1:n
        if ~visited(start)
            % Fix the starting triangle as-is and propagate constraints
            queue = start;
            visited(start) = true;

            while ~isempty(queue)
                t = queue(1);
                queue(1) = [];

                tri = facet_nds_out(t,:);
                edges_t = [tri([1 2]); tri([2 3]); tri([3 1])]; % oriented edges of t

                % For each edge of t, enforce opposite orientation in neighbor
                for k = 1:3
                    ut = edges_t(k,1); vt = edges_t(k,2);           % t's oriented edge ut->vt
                    key = edge_key(ut, vt);
                    if ~isKey(edgeMap, key)
                        continue; % should not happen
                    end
                    tris_sharing = edgeMap(key);
                    % There can be one (boundary) or two triangles for this edge
                    % Process the neighbor (if present)
                    for idx = 1:numel(tris_sharing)
                        j = tris_sharing(idx);
                        if j == t, continue; end

                        if ~visited(j)
                            % Reorder j so that it has the edge oriented opposite: vt->ut
                            nodes_j = facet_nds_out(j,:);
                            % Identify the third node r in triangle j
                            mask_uv = ismember(nodes_j, [ut vt]);
                            if nnz(mask_uv) ~= 2
                                % This should not happen if input is consistent
                                warning('Triangle %d does not contain the expected shared edge (%d,%d).', j, ut, vt);
                                continue;
                            end
                            r = nodes_j(~mask_uv); % the remaining node

                            % Set the order explicitly to enforce vt->ut as a consecutive edge
                            facet_nds_out(j,:) = [vt, ut, r];

                            visited(j) = true;
                            queue(end+1) = j;
                        else
                            % Already oriented: verify it has vt->ut as a consecutive edge
                            nodes_j = facet_nds_out(j,:);
                            if ~has_oriented_edge(nodes_j, vt, ut)
                                % Conflict detected: changing j now could break already-set neighbors,
                                % so we report and leave as-is.
                                warning(['Orientation conflict on shared edge {%d,%d} between triangles %d and %d. ', ...
                                         'Triangle %d does not have the required reversed edge.'], ...
                                         min(ut,vt), max(ut,vt), t, j, j);
                            end
                        end
                    end
                end
            end
        end
    end
%list of all edges and associated faces
eds = [facet_nds_out(:,1), facet_nds_out(:,2)
    facet_nds_out(:,2), facet_nds_out(:,3)
    facet_nds_out(:,3), facet_nds_out(:,1)];
ed_fcs = repmat((1:size(facet_nds_out, 1))', [3, 1]);
end

% -------------------- Helper functions --------------------
function key = edge_key(u, v)
    % Build undirected edge key "min-max"
    a = min(u, v);
    b = max(u, v);
    key = sprintf('%d-%d', a, b);
end

function tf = has_oriented_edge(row, s, t)
    % True if row has consecutive oriented edge s->t in positions (1,2), (2,3), or (3,1)
    tf = (row(1)==s && row(2)==t) || ...
         (row(2)==s && row(3)==t) || ...
         (row(3)==s && row(1)==t);
end


% %gets matrix of facet nodes into consistent ordering (say clockwise) around
% %each facet, so surface normals of adjacent facets are in same direction
% 
% if size(facet_nds, 2) ~= 3
%     error('Facet nodes must have 3 columns (triangular facets')
% end
% 
% 
% %Process all the facet edges
% processed = zeros(size(facet_nds, 1), 1);
% flipped_out = zeros(size(facet_nds, 1), 1);
% matched = zeros(size(facet_nds));
% i = 1;
% while 1
%     processed(i) = 1;
%     %loop edges of current facet
%     tmp_ed1 = [facet_nds(i,[1,2]);facet_nds(i,[2, 3]); facet_nds(i,[3,1])];
%     tmp_ed1 = fn_eds_for_facet(facet_nds(i, :));
%     % tmp_ed1 = eds(ed_fcs == i, :);
%     for j = 1:3
%         %find other facet containing current edge
%         k = find(any(ismember(facet_nds, tmp_ed1(j, 1)), 2) & any(ismember(facet_nds, tmp_ed1(j, 2)), 2) & ~processed);
%         if ~isempty(k)
%             %get edges of that facet
%             tmp_ed2 = [facet_nds(k,[1,2]);facet_nds(k,[2, 3]); facet_nds(k,[3,1])];
%             tmp_ed2 = fn_eds_for_facet(facet_nds(k, :));
%             processed(k) = 1;
%             matched(i, j) = 1;
%             %see if current edge is same order in second facet and flip current edge if
%             %nesc (they need to be in opposite orders!!!!)
%             if any(ismember(tmp_ed2, tmp_ed1(j, :), 'rows'))
%                 facet_nds(k,:) = fliplr(facet_nds(k,:));
%                 flipped_out(k) = 1;
% 
%             end
%         end
%     end
%     i = find(~processed);
%     if isempty(i)
%         %all facets processed
%         break
%     else
%         %move to next unprocessed facet
%         i = i(1);
%     end
% end
% 
% %list of all edges and associated faces
% eds = [facet_nds(:,1), facet_nds(:,2)
%     facet_nds(:,2), facet_nds(:,3)
%     facet_nds(:,3), facet_nds(:,1)];
% ed_fcs = repmat((1:size(facet_nds, 1))', [3, 1]);
% 
% 
% 
% 
% 
% % [~,ia,ic] = unique(sort(eds, 2), 'rows');
% % occs = accumarray(ic,1);
% % if any(occs) > 2
% %     error('More than two facets at same edge')
% % end
% % occs = find(occs > 1);
% % for i = 1:numel(occs)
% %     j = find(ic == occs(i));
% %     if all(eds(j(1), :) == eds(j(2), :))
% %         %Flip order of nodes in associated face
% %         flip_face = ed_fcs(j(1));
% %         facet_nds(flip_face,:) = fliplr(facet_nds(flip_face,:));
% %         eds(j(1), :) = fliplr(eds(j(1), :));
% %     end
% % end
% end
% 
% function eds = fn_eds_for_facet(facet_nds)
% eds = [facet_nds(:, [1, 2]); facet_nds(:, [2, 3]); facet_nds(:, [3, 1])];
% end