function [max_velocity, min_velocity] = fn_estimate_max_min_vels(matls)
%USAGE
%   [max_velocity, min_velocity] = fn_estimate_max_min_vels(matls)
%AUTHOR
%   Paul Wilcox (2025)
%SUMMARY
%   Estimates minimum and maximum speeds of sound possible given cell array
%   of materials
%INPUTS
%   matls - cell array of materials with fields D (stiffness matrix
%   [solids] or bulk modulus [fluids]), rho (density)
%OUTPUTS
%   max_velocity, min_velocity - as per variable names
%NOTES
%   This does a quick and easy calculation that is correct for isotropic
%   solids and fluids (in that it returns longitudinal and shear velocities
%   for former, and pressure velocity for latter). For anisotropic
%   materials it will return plausible values but they will not correspond
%   to specific points on slowness surfaces.
%--------------------------------------------------------------------------

%Deal with legacy v2 code where materials is structure array rather than cell array
if isstruct(matls)
    matls = arrayfun(@(x) x, matls, 'UniformOutput', false);
end

v = [0,0];
j = 1;
for i = 1:numel(matls)
    if isempty(matls{i})
        continue
    end
    s = [];
    if isfield(matls{i}, 'stiffness_matrix')
        s = matls{i}.stiffness_matrix(:);
    end
    if isfield(matls{i}, 'D')
        s = matls{i}.D(:);
    end
    if isfield(matls{i}, 'density') %deal with legacy naming
        rho = matls{i}.density;
    else
        rho = matls{i}.rho;
    end

    if ~isempty(s)
        s = s(abs(s) > 0);
        v(j,1) = sqrt(min(s) / rho);
        v(j,2) = sqrt(max(s) / rho);
        j = j + 1;
        % s = s(abs(s) > 0);
        % v(j,1) = sqrt(min(s) / rho);
        % v(j,2) = sqrt(max(s) / rho);
        % j = j + 1;
    end
end
min_velocity = min(v(:,1));
max_velocity = max(v(:,2));
end