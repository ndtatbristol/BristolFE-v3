function [mod, matls, el_types, steps, fe_options, params] = mod_2d_basic(params)
%This is an example of a parametric model description function. The model
%is built based on the parameters supplied. All the parameters have default
%values at the top of the file, so any that are not specified take the
%default values. This provides a clean way of producing re-usable 
%parameteric models and leads to tidy scripts to actually run them 
%(possibly within the context of a parametric study). 
%This particular model is used in the following scripts:
% - ex_2d_basic_using_mod.m
% - ex_2d_bristolFE_pogo_comparison.m

default_params.model_size = 20e-3;
default_params.els_per_wavelength = 10;

%Source details (expressed here in terms of fractions of model_size)
default_params.source_position_as_fractions = [1 / 5, 1 / 7];
default_params.source_direction = 2; %2 means the second DoF. Possible DoF for solids are 1: x, 2: y, 3: z

%Monitoring details (expressed here in terms of fractions of model_size)
default_params.monitor_position_as_fractions = [1 / 3, 1 / 4];
default_params.monitor_direction = 2; %2 means the second DoF. Possible DoF for solids are 1: x, 2: y, 3: z

%Material properties
default_params.matl_name = 'aluminium';

%Details of input signal applied at source
default_params.centre_freq = 5e6;
default_params.no_cycles = 4;
%Run for long enough for longitudinal waves to travel this many lengths of model
default_params.max_time_multiplier = 3; 

%Element shape to use (tri or quad since this is 2d model)
default_params.element_shape = 'tri'; 

%Solver options - specify how ofter field output is produced to use in
%animation
default_params.fe_options.field_output_every_n_frames = 20;

%--------------------------------------------------------------------------
params.fe_options = fn_set_default_fields(params.fe_options, default_params.fe_options);
params = fn_set_default_fields(params, default_params);
fe_options = params.fe_options;
el_types = fn_2d_el_types();

switch params.element_shape
    case {'tri', 'triangular'}
        el_typ_to_use_for_solid = 'CPE3';
        el_typ_to_use_for_fluid = 'AC2D3';
    case {'quad', 'rect', 'quadrilateral', 'rectangle', 'rectangular'}
        el_typ_to_use_for_solid = 'CPE4';
        el_typ_to_use_for_fluid = 'AC2D4';
    otherwise
        error('Unknown element shape')
end


%BUILD THE MODEL USING THE PARAMETERS GIVEN ABOVE
%Polygon vertices define shape of model (in this case a square with side
%length = model_size)
bdry_pts = [
    0,          0 
    params.model_size, 0 
    params.model_size, params.model_size
    0,          params.model_size];

%Define the material
%A cell array with an entry for each material used in the model is required.
%In this case there is only one material, but the index of this material in 
%the cell array is given a name, matl_i, so you can see where it is used
%later when elements in the model are assigned to elements
matl_i = 1; %material index is given a name so you can see where it appears later
matls{matl_i} = fn_material_library(params.matl_name);

%Work out element size and time step
el_size = fn_get_suitable_el_size(matls, params.centre_freq, params.els_per_wavelength);
time_step = fn_get_suitable_time_step(matls, el_size);

%Associate each element with a material index and element type index
switch params.element_shape
    case {'tri', 'triangular'}
        el_typ_to_use_for_solid = 'CPE3';
    case {'quad', 'rect', 'quadrilateral', 'rectangle', 'rectangular'}
        el_typ_to_use_for_solid = 'CPE4';
    otherwise
        error('Unknown element shape')
end

%Create the nodes and elements of the mesh
mod = fn_2d_structured_mesh(bdry_pts, el_size, el_typ_to_use_for_solid);

mod.el_mat_i(:) = matl_i;
mod.el_typ_i(:) = find(strcmp(el_types, el_typ_to_use_for_solid)); %extracts the index of the chosen element type from the cell array of possible element types

%Identify node closest to desired source location
source_position = params.model_size * params.source_position_as_fractions;
steps{1}.load.frc_nds = fn_find_node_nearest_to_point(mod.nds, source_position, el_size);
steps{1}.load.frc_dfs = ones(size(steps{1}.load.frc_nds)) * params.source_direction;

%Provide the time signal for the loading
[vel, ~] = fn_estimate_max_min_vels(matls{matl_i});
max_time = params.model_size / vel * params.max_time_multiplier;
steps{1}.load.time = 0: time_step: max_time;
steps{1}.load.frcs = fn_gaussian_pulse(steps{1}.load.time, params.centre_freq, params.no_cycles);

%Say where the displacement should be monitored
monitor_position = params.model_size * params.monitor_position_as_fractions;
steps{1}.mon.nds = fn_find_node_nearest_to_point(mod.nds, monitor_position, el_size);
steps{1}.mon.dfs =  ones(size(steps{1}.mon.nds)) * params.monitor_direction;

end