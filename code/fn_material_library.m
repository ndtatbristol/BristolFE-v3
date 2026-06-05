function matl = fn_material_library(name)
%USAGE
%   matl = fn_material_library(name)
%SUMMARY
%   Returns material structure for material with given name if found in
%   libary
%INPUTS
%   name - name of material to look for
%RETURNS
%   matls - structure in correct form for usin BristolFE
%--------------------------------------------------------------------------

switch lower(name)
    case {'aluminium', 'al', 'aluminum'}
        longitudinal_velocity = 6300;
        shear_velocity = 3150;
        density = 2700;
        matl = fn_matl_isotropic_solid_defined_by_velocities(name, longitudinal_velocity, shear_velocity, density);
    case {'steel', 'st'}
        longitudinal_velocity = 5940;
        shear_velocity = 3220;
        density = 7850;
        matl = fn_matl_isotropic_solid_defined_by_velocities(name, longitudinal_velocity, shear_velocity, density);
    case {'perspex', 'plexiglass'}
        longitudinal_velocity = 2730;
        shear_velocity = 1345;
        density = 1190;
        matl = fn_matl_isotropic_solid_defined_by_velocities(name, longitudinal_velocity, shear_velocity, density);
    case {'water', 'h2o'}
        velocity = 1500;
        density = 1000;
        matl = fn_matl_fluid_defined_by_velocity(name, velocity, density);
    otherwise
        error(['Material "', name, '" not found in material library']);
end
