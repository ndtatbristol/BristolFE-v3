function matl = fn_matl_isotropic_solid_defined_by_stiffness(name, youngs_modulus, poissons_ratio, density, varargin)
%USAGE
%   matl = fn_matl_isotropic_solid_defined_by_stiffness(name, youngs_modulus, poissons_ratio, density [, col, options])
%SUMMARY
%   Returns material structure based on specifed youngs_modulus, 
%   poissons_ratio, and density.
%INPUTS
%   name - name of material
%   youngs_modulus, poissons_ratio, density - as per variable names
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

default_options.min_colour = [0.75, 0.75, 0.75];
default_options.max_colour = [0.1, 0.1, 0.1];
default_options.zmin = 2e6;
default_options.zmax = 50e6;
options = fn_set_default_fields(options, default_options);

% base_colour_hsv = rgb2hsv(options.base_colour);

matl.rho = density;
matl.D = fn_isotropic_stiffness_matrix(youngs_modulus, poissons_ratio); 
matl.name = name;

z = sqrt(mean(matl.D, 'all') * matl.rho); %sort of proxy for acoustic impedance, used to generate default colour and name
if ~isempty(col)
    matl.col = col;
else
    %Col
    zc = max(z, options.zmin);
    zc = min(zc, options.zmax);
    % matl.col = hsv2rgb([base_colour_hsv(1), base_colour_hsv(2), (options.zmax - zc) / (options.zmax - options.zmin)]);
    matl.col = interp1([0,1], [options.min_colour; options.max_colour], (zc - options.zmin) / (options.zmax - options.zmin), 'linear');
end
