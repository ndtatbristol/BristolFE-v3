function matl = fn_matl_fluid_defined_by_velocity(name, velocity, density, varargin)
%USAGE
%   matl = fn_matl_fluid_defined_by_velocity(name, velocity, density [, col, options])
%SUMMARY
%   Returns material structure for fluid based on specifed velocity and density.
%INPUTS
%   name - name of material
%   velocity, density - as per variable names
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
bulk_modulus = velocity ^ 2 * density;
matl = fn_matl_fluid_defined_by_bulk_modulus(name, bulk_modulus, density, col, options);

end
