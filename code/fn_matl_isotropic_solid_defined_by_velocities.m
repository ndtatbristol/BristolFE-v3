function matl = fn_matl_isotropic_solid_defined_by_velocities(name, longitudinal_velocity, shear_velocity, density, varargin)
%USAGE
%   matl = fn_matl_isotropic_solid_defined_by_velocities(name, longitudinal_velocity, shear_velocity, density, [, col, options])
%SUMMARY
%   Returns material structure based on specifed velocities and density.
%INPUTS
%   name - name of material
%   longitudinal_velocity, shear_velocity, density - as per variable names
%   [col - RGB colour of material to use in plots. Will be generated
%   automatically if not specified]
%   [options - options structure that can be used to specify base colour,
%   zmin and zmax values for calculating actual colour
%OUTPUT
%   matl structure with fields D, rho, name, col
%--------------------------------------------------------------------------
if numel(varargin) >= 1
    col = varargin{1};
else
    col = [];
end
if numel(varargin) >= 2
    options = varargin{2};
else
    options = [];
end

[youngs_modulus, poissons_ratio] = fn_stiffness_from_velocities_and_density(longitudinal_velocity, shear_velocity, density);

matl = fn_matl_isotropic_solid_defined_by_stiffness(name, youngs_modulus, poissons_ratio, density, col, options);
end