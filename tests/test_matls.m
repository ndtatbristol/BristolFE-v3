clear
close all

matls{1} = fn_matl_isotropic_solid_defined_by_stiffness('Steel', 210e9, 0.3, 8900);
matls{2} = fn_matl_isotropic_solid_defined_by_stiffness('Aluminium', 70e9, 1/3, 2700);
matls{3} = fn_matl_isotropic_solid_defined_by_velocities('St', 5900, 3150, 8900);
matls{4} = fn_matl_isotropic_solid_defined_by_velocities('Al', 6300, 3150, 2700);
matls{5} = fn_matl_fluid_defined_by_velocity('water', 1500, 1000);
matls{6} = fn_matl_fluid_defined_by_bulk_modulus('water', 2.2500e+09, 1000);

for i = 1:numel(matls)
    [max_velocity, min_velocity] = fn_estimate_max_min_vels(matls{i});
    fprintf('\nMaterial %i\n', i);
    fprintf(['Name:', matls{i}.name, '\n']);
    fprintf('Max vel: %.1f m/s\n', max_velocity);
    fprintf('Min vel: %.1f m/s\n', min_velocity);
end

[max_velocity, min_velocity] = fn_estimate_max_min_vels(matls);
fprintf('\n\nAll materials\n');
fprintf('Max vel: %.1f m/s\n', max_velocity);
fprintf('Min vel: %.1f m/s\n', min_velocity);

centre_frequency = 5e6;
[max_wavelength, min_wavelength] = fn_estimate_max_min_wavelengths(matls, centre_frequency);
fprintf('Max wavelength: %.3f mm\n', max_wavelength * 1e3);
fprintf('Min wavelength: %.3f mm\n', min_wavelength * 1e3);

