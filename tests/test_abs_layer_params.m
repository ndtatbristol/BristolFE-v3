clear all;
close all;

%Material properties
matl_longitudinal_velocity = 1;
matl_shear_velocity = 0.5;
matl_density = 1;
matl_name = 'magic';

centre_freq = 1;

lambda = matl_longitudinal_velocity / centre_freq;

%should check it works over large range of scales
model_size = 10 * lambda;
abs_bdry_thicknesses = [1:0.5:3] * lambda;
max_damping = centre_freq * linspace(8, 40, 20);

% abs_bdry_thicknesses = [1, 3] * lambda;
%Elements per wavelength which is used to determine element size (more
%elements per wavelength = more accurate model but higher computational cost)
els_per_wavelength = 10;

%Source details (expressed here in terms of model_size)
source_position = [model_size * 0.45, model_size * 0.45];
source_direction = 2; %2 means the second DoF. Possible DoF for solids are 1: x, 2: y, 3: z

%Monitoring details (expressed here in terms of model_size)
monitor_position = [model_size * 0.55, model_size * 0.55];
monitor_direction = 2; %2 means the second DoF. Possible DoF for solids are 1: x, 2: y, 3: z


%Details of input signal applied at source
no_cycles = 4;
%Run for long enough for longitudinal waves to travel 3 lengths of model
max_time = 3 * model_size / matl_longitudinal_velocity;

%Say which element type will be used (currently there is only one choice for
%a solid material in 2D model, but in the future there may be more options to choose
%from, e.g. quadrilateral elements, second order elements)
el_typ_to_use_for_solid = 'CPE3';

%Solver options - specify how ofter field output is produced to use in
%animation
fe_options.field_output_every_n_frames = 5;

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
mod = fn_2d_structured_mesh_triangular_els(bdry_pts, el_size);

%A cell array that includes all element types used in a model is required
%(same idea as the cell array of materials that contains all materials used
%in the model). The function below produces a list of all available 2d
%elements, which is fine (it doesn't matter that some won't be used)
el_types = fn_2d_el_types();

%Associate each element with a material index and element type index
mod.el_mat_i(:) = matl_i;
mod.el_typ_i(:) = find(strcmp(el_types, el_typ_to_use_for_solid)); %extracts the index of the chosen element type from the cell array of possible element types

%Identify node closest to desired source location
steps{1}.load.frc_nds = fn_find_node_nearest_to_point(mod.nds, source_position, el_size);
steps{1}.load.frc_dfs = ones(size(steps{1}.load.frc_nds)) * source_direction;

%Provide the time signal for the loading
steps{1}.load.time = 0: time_step:  max_time;
steps{1}.load.frcs = fn_gaussian_pulse(steps{1}.load.time, centre_freq, no_cycles);

%Say where the displacement should be monitored
steps{1}.mon.nds = fn_find_node_nearest_to_point(mod.nds, monitor_position, el_size);
steps{1}.mon.dfs =  ones(size(steps{1}.mon.nds)) * monitor_direction;


for i = 1:numel(abs_bdry_thicknesses)
    abs_bdry_thickness = abs_bdry_thicknesses(i);
    abs_bdry_pts = [
        abs_bdry_thickness,                 abs_bdry_thickness
        model_size - abs_bdry_thickness,  abs_bdry_thickness
        model_size - abs_bdry_thickness,  model_size - abs_bdry_thickness
        abs_bdry_thickness,                 model_size - abs_bdry_thickness];

    mod = fn_2d_add_absorbing_layer(mod, abs_bdry_pts, abs_bdry_thickness);

    %--------------------------------------------------------------------------
    %RUN THE MODEL
    % max_damping = centre_freq * linspace(5, 10, 10);
    %fe_options.max_stiffness_reduction = 0.1;
    for j = 1:numel(max_damping)
        fe_options.max_damping = max_damping(j);
        res = fn_FE_entry_point(mod, matls, el_types, steps, fe_options);
        tmp = sum(abs(res{1}.fld));
        if (i == 1) && (j == 1)
            KE = zeros(numel(tmp), numel(max_damping), numel(abs_bdry_thicknesses));
        end
        KE(:, j, i) = tmp;
    end
end
KE = KE / max(KE, [], 'all');

figure
for i = 1:numel(abs_bdry_thicknesses)
    hold on
    plot(max_damping, 10*log10(KE(end,:, i)));
end
return

%--------------------------------------------------------------------------
%SHOW THE RESULTS

%Plot history output at monitoring node
figure;
semilogy(steps{1}.load.time, abs(res{1}.dsps) / max(abs(res{1}.dsps)));
ylim([0.0001, 1])
xlabel('Time (s)')
title(sprintf('Max damp: %.2e, thicknes: %.2e mm', fe_options.max_damping, abs_bdry_thickness));
return
%Animate field output result
figure;
%Exactly the same function as that used to show the geometry is used to
%produce the plot for animation. It returns a handle to the patches representing
%elements and it is the colours of these that are animated by fn_run_animation
display_options = [];
h_patch = fn_show_geometry(mod, matls, display_options);
anim_options.db_range = [-40, 0];
fn_run_animation(h_patch, res{1}.fld, anim_options);
