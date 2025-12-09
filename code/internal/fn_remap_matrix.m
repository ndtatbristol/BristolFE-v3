function m2 = fn_remap_matrix(m1, map, varargin)
%REMAP_MATRIX Remap integer values in m1 using a 2-col mapping [old new].
%   m2 = fn_remap_matrix(m1, map)
%   m2 = fn_remap_matrix(m1, map, 'MissingPolicy', policy, 'MissingValue', mv)
%
% Inputs:
%   m1  : matrix of integers (any shape)
%   map : Kx2 integer matrix; map(:,1)=old, map(:,2)=new
%OR
%   map : Kx1 integer matrix of new where old is assumed to be 1:K
%
% Name-Value:
%   'MissingPolicy' : 'keep' | 'nan' (default) | 'error' | 'fill'
%   'MissingValue'  : value used when MissingPolicy='fill' (default=NaN)
%
% Output:
%   m2 : matrix with values replaced per map
%
% Notes:
%   - Uses a vectorized ismember-based method (works for any integer keys).
%   - For speed with dense positive keys, you can swap to a LUT approach.

    p = inputParser;
    addParameter(p, 'MissingPolicy', 'nan', @(s) any(strcmpi(s, {'keep','nan','error','fill'})));
    addParameter(p, 'MissingValue', NaN);
    parse(p, varargin{:});
    policy = lower(p.Results.MissingPolicy);
    mv = p.Results.MissingValue;

    % Basic checks
    if size(map,2) == 1
        map = [(1:numel(map))', map];
    end
    if size(map,2) ~= 2
        error('map must be Kx2.');
    end
    if ~isnumeric(m1) || ~isnumeric(map)
        error('m1 and map must be numeric.');
    end

    [tf, loc] = ismember(m1, map(:,1));
    switch policy
        case 'keep'
            m2 = m1;
            m2(tf) = map(loc(tf), 2);

        case 'nan'
            m2 = nan(size(m1));
            m2(tf) = map(loc(tf), 2);

        case 'fill'
            m2 = ones(size(m1), 'like', m1) * mv;
            % ensure the right class if mv is double and m1 is integer
            if ~isa(m2, class(map(:,2)))
                m2 = cast(m2, class(map(:,2)));
            end
            m2(tf) = map(loc(tf), 2);

        case 'error'
            if any(~tf(:))
                missing = unique(m1(~tf));
                error('Missing keys in map for values: %s', mat2str(missing(:)'));
            end
            m2 = map(loc, 2);

        otherwise
            error('Unknown MissingPolicy: %s', policy);
    end
end
