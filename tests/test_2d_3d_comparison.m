clear
close all
%This builds a 2D model and a 3D model one element thick with mirror BCs -
%these should give same results

model_size = 20e-3;

%Elements per wavelength which is used to determine element size (more 
%elements per wavelength = more accurate model but higher computational cost)
els_per_wavelength = 10;

%Source details (expressed here in terms of model_size)
source_position_2d = [model_size / 5, model_size / 7];
source_direction = 2; %2 means the second DoF. Possible DoF for solids are 1: x, 2: y, 3: z

%Monitoring details (expressed here in terms of model_size)
monitor_position_2d = [model_size / 3, model_size / 4];
monitor_direction = 2; %2 means the second DoF. Possible DoF for solids are 1: x, 2: y, 3: z

%Material properties
matl_longitudinal_velocity = 6300;
matl_shear_velocity = 3150;
matl_density = 2700;
matl_name = 'aluminium';

%Details of input signal applied at source
centre_freq = 5e6;
no_cycles = 4;
%Run for long enough for longitudinal waves to travel 3 lengths of model
max_time = 3 * model_size / matl_longitudinal_velocity; 

%Say which element type will be used (currently there is only one choice for
%a solid material in 2D model, but in the future there may be more options to choose
%from, e.g. quadrilateral elements, second order elements)
el_typ_to_use_for_solid_2d = 'CPE4'; 
el_typ_to_use_for_solid_3d = 'C3D8'; 

%Solver options - specify how ofter field output is produced to use in
%animation
fe_options.field_output_every_n_frames = 5;
show_geom = 0;
%--------------------------------------------------------------------------
%THE ACTUAL CODE STARTS HERE
%--------------------------------------------------------------------------

%SET UP THE MODEL
%Add path to BristolFE functions in case not already on path
addpath(genpath('../code')); 

%BUILD THE MODEL USING THE PARAMETERS GIVEN ABOVE
%Polygon vertices define shape of model (in this case a square with side
%length = model_size)
bdry_pts = [
    0,          0 
    model_size, 0 
    model_size, model_size
    0,          model_size];

%Define the material
%A cell array with an entry for each material used in the model is required.
%In this case there is only one material, but the index of this material in 
%the cell array is given a name, matl_i, so you can see where it is used
%later when elements in the model are assigned to elements
matl_i = 1; %material index is given a name so you can see where it appears later
matls{matl_i} = fn_matl_isotropic_solid_defined_by_velocities(matl_name, matl_longitudinal_velocity, matl_shear_velocity, matl_density);

%Work out element size and time step
el_size = fn_get_suitable_el_size(matls, centre_freq, els_per_wavelength);
time_step = fn_get_suitable_time_step(matls, el_size);

%Create the nodes and elements of the mesh
mod_2d = fn_2d_structured_mesh_rectangular_els(bdry_pts, el_size);
crnr_pts = [0, 0, 0
    bdry_pts(3, :), el_size];
mod_3d = fn_3d_structured_mesh_hexahedral_els(crnr_pts, el_size);

%A cell array that includes all element types used in a model is required 
%(same idea as the cell array of materials that contains all materials used 
%in the model). The function below produces a list of all available 2d
%elements, which is fine (it doesn't matter that some won't be used)
el_types_2d = fn_2d_el_types(); 
el_types_3d = fn_3d_el_types(); 

%Associate each element with a material index and element type index
mod_2d.el_mat_i(:) = matl_i;
mod_2d.el_typ_i(:) = find(strcmp(el_types_2d, el_typ_to_use_for_solid_2d)); %extracts the index of the chosen element type from the cell array of possible element types

mod_3d.el_mat_i(:) = matl_i;
mod_3d.el_typ_i(:) = find(strcmp(el_types_3d, el_typ_to_use_for_solid_3d)); %extracts the index of the chosen element type from the cell array of possible element types

%Identify node closest to desired source location
steps_2d{1}.load.frc_nds = fn_find_node_nearest_to_point(mod_2d.nds, source_position_2d, el_size);
steps_2d{1}.load.frc_dfs = ones(size(steps_2d{1}.load.frc_nds)) * source_direction;

steps_3d{1}.load.frc_nds = [...
    fn_find_node_nearest_to_point(mod_3d.nds, [source_position_2d, 0], el_size), ...
    fn_find_node_nearest_to_point(mod_3d.nds, [source_position_2d, el_size], el_size)];
steps_3d{1}.load.frc_dfs = ones(size(steps_3d{1}.load.frc_nds)) * source_direction;




%Provide the time signal for the loading
steps_2d{1}.load.time = 0: time_step: max_time;
steps_2d{1}.load.frcs = fn_gaussian_pulse(steps_2d{1}.load.time, centre_freq, no_cycles);
steps_3d{1}.load.time = 0: time_step: max_time;
steps_3d{1}.load.frcs = fn_gaussian_pulse(steps_3d{1}.load.time, centre_freq, no_cycles) * el_size / 2;
%all nodes need to have zero disp in z
steps_3d{1}.load.dsp_nds = 1:size(mod_3d.nds, 1);
steps_3d{1}.load.dsp_dfs = ones(size(steps_3d{1}.load.dsp_nds)) * 3;
steps_3d{1}.load.dsps = zeros(size(steps_3d{1}.load.time));

%Say where the displacement should be monitored
steps_2d{1}.mon.nds = fn_find_node_nearest_to_point(mod_3d.nds(:,1:2), monitor_position_2d, el_size);
steps_2d{1}.mon.dfs =  ones(size(steps_2d{1}.mon.nds)) * monitor_direction;
steps_3d{1}.mon.nds = [...
    fn_find_node_nearest_to_point(mod_3d.nds, [monitor_position_2d, 0], el_size), ...
    fn_find_node_nearest_to_point(mod_3d.nds, [monitor_position_2d, el_size], el_size)];
steps_3d{1}.mon.dfs = ones(size(steps_3d{1}.mon.nds)) * monitor_direction;

%At this point the model is completely defined by mod, matls, el_types, and 
%steps. For example
%
%       mod.nds: [128205×2 double]
%      mod.els: [254976×3 double]
% mod.el_mat_i: [254976×1 double]
% mod.el_typ_i: [254976×1 double]
% mod.el_abs_i: [254976×1 double]
% 
% matls{1}.rho: 8900
%   matls{1}.D: [6×6 double]
%matls{1}.name: 'steel'
% matls{1}.col: [0.2040 0.2040 0.2040]
%
%     el_types: {'CPE3', 'AC2D3', 'ASI2D2'}
%
% steps{1}.load: [1×1 struct]
% steps{1}.mon: [1×1 struct]

%--------------------------------------------------------------------------
%SHOW THE MESH
%This is an optional step, but sometimes you might want to look at the mesh 
%before running the model to make sure it looks correct.
%The display_options parameter includes the option for sets o
%of nodes to be marked. Here it is used to show the source and monitor
%nodes.
if show_geom
    figure;
    subplot(1,2,1);
    display_options.node_sets_to_plot(1).nd = steps_2d{1}.load.frc_nds;
    display_options.node_sets_to_plot(1).col = 'r.';
    display_options.node_sets_to_plot(2).nd = steps_2d{1}.mon.nds;
    display_options.node_sets_to_plot(2).col = 'g.';
    fn_show_geometry(mod_2d, matls, el_types_2d, display_options);
    subplot(1,2,2);
    display_options.node_sets_to_plot(1).nd = steps_3d{1}.load.frc_nds;
    display_options.node_sets_to_plot(1).col = 'r.';
    display_options.node_sets_to_plot(2).nd = steps_3d{1}.mon.nds;
    display_options.node_sets_to_plot(2).col = 'g.';
    fn_show_geometry(mod_3d, matls, el_types_3d, display_options);
    drawnow
end
%--------------------------------------------------------------------------
%RUN THE MODEL

%This is where the model actually gets executed.
res_2d = fn_FE_entry_point(mod_2d, matls, el_types_2d, steps_2d, fe_options);
res_3d = fn_FE_entry_point(mod_3d, matls, el_types_3d, steps_3d, fe_options);

%--------------------------------------------------------------------------
%SHOW THE RESULTS

%Plot history output at monitoring node
figure;
plot(steps_2d{1}.load.time, res_2d{1}.dsps, 'b');
hold on;
plot(steps_3d{1}.load.time, res_3d{1}.dsps, 'r:');
plot(steps_3d{1}.load.time, sum(res_3d{1}.dsps, 1), 'r');
xlabel('Time (s)')
legend('2D', '3D')
return
%Animate field output result
figure;
%Exactly the same function as that used to show the geometry is used to
%produce the plot for animation. It returns a handle to the patches representing
%elements and it is the colours of these that are animated by fn_run_animation
h_patch = fn_show_geometry(mod_2d, matls, el_types_2d, display_options);
anim_options.fld_time = res_2d{1}.fld_time;
fn_run_animation(h_patch, res_2d{1}.fld, anim_options);
