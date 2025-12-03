function el_size = fn_get_suitable_el_size(matls, nominal_cent_freq, els_per_wavelength)
%USAGE
%   el_size = fn_get_suitable_el_size(matls, nominal_cent_freq, els_per_wavelength)
%AUTHOR
%   Paul Wilcox (2025)
%SUMMARY
%    Estimated element size to achieve the desired number of elements per
%    wavelength for the slowest wave possible in any of the materials
%    given. For isotropic solids, the slowest wave will be the shear mode,
%    so for solid models, this function returns an element size based on 
%    the number of elements per SHEAR wavelength at the centre frequency.
%INPUT
%   matls - cell array of materials
%   nominal_cent_freq - centre frequency at which to compute wavelength
%   els_per_wavelength - minimum number of elements required per wavelength
%OUTPUT
%   el_size - the calculated element size
%--------------------------------------------------------------------------

[~, min_vel] = fn_estimate_max_min_vels(matls);
lambda_min = min_vel / nominal_cent_freq;

el_size = lambda_min / els_per_wavelength;
end
