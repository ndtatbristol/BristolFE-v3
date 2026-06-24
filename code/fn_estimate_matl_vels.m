function [max_velocity, min_velocity] = fn_estimate_matl_vels(matl)
%USAGE
%   [max_velocity, min_velocity] = fn_estimate_matl_vels(matls)
%AUTHOR
%   Paul Wilcox (2025)
%SUMMARY
%   Estimates minimum and maximum speeds of sound possible for given
%   material. For isotropic elastic solids these correspond to shear and
%   longitudinal wave speeds. For fluids they will be same and correspond
%   to bulk pressure wave speed. For anisotropic solids, who knows?
%INPUTS
%   matl - material structure with fields D (stiffness matrix
%   [solids] or bulk modulus [fluids]), rho (density)
%OUTPUTS
%   max_velocity, min_velocity - as per variable names
%NOTES
%   This does a quick and easy calculation that is correct for isotropic
%   solids and fluids (in that it returns longitudinal and shear velocities
%   for former, and pressure wave velocity for latter). For anisotropic
%   materials it will return plausible values but they will not correspond
%   to specific points on slowness surfaces.
%--------------------------------------------------------------------------

if isfield(matl, 'stiffness_matrix')
    s = matl.stiffness_matrix(:);
elseif isfield(matl, 'D')
    s = matl.D(:);
else
    error('No stiffness field found in material')
end
if isfield(matl, 'density') %deal with legacy naming
    rho = matl.density;
elseif isfield(matl, 'rho')
    rho = matl.rho;
else
    error('No density field found in material')
end

s = s(abs(s) > 0);
min_velocity = sqrt(min(s) / rho);
max_velocity = sqrt(max(s) / rho);
end