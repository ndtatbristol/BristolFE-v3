function [max_wavelength, min_wavelength] = fn_estimate_max_min_wavelengths(matls, centre_frequency)
%USAGE
%   [max_wavelength, min_wavelength] = fn_estimate_max_min_wavelengths(matls, centre_frequency)
%AUTHOR
%   Paul Wilcox (2025)
%SUMMARY
%   Estimates minimum and maximum wavelengths possible given cell array
%   of materials
%INPUTS
%   matls - cell array of materials with fields D (stiffness matrix
%   [solids] or bulk modulus [fluids]), rho (density)
%OUTPUTS
%   max_wavelength, min_wavelength - as per variable names
%NOTES
%   This does a quick and easy calculation that is correct for isotropic
%   solids and fluids (in that it returns longitudinal and shear
%   wavelengths for former, and pressure wavelength for latter). For anisotropic
%   materials it will return plausible values but they will not correspond
%   to specific points on slowness surfaces.
%--------------------------------------------------------------------------
[max_velocity, min_velocity] = fn_estimate_max_min_vels(matls);
tmp = [max_velocity, min_velocity] / centre_frequency;
max_wavelength = tmp(1);
min_wavelength = tmp(2);
end