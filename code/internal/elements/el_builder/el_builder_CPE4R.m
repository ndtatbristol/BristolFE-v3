%test new element formulation method
clear all
restoredefaultpath;
addpath(genpath('..\..\..\..\code'));

%CPE4
% nds_in_nat_coords = [
%     -1, -1 
%      1, -1 
%      1,  1 
%     -1,  1];
nds_in_nat_coords = [
    -1, -1 
    -1,  1 
     1,  1 
     1, -1];
sf_powers = [
    0, 0 
    1, 0 
    0, 1 
    1, 1];

% reduced integration
gauss_pts = [
    0, 0];
gauss_weights = 4; 

no_dfs = 3;
solid_or_fluid = 'solid';
new_el_type = 'CPE4R';
ref_el_type = 'CPE4R';
no_dims = 2;

[K, M, detJ, loc_nd, loc_df] = fn_symbolic_K_and_M_matrices(nds_in_nat_coords, gauss_pts, gauss_weights, sf_powers, no_dfs, solid_or_fluid);

C = sym(zeros(size(K)));

fn_create_element_matrix_file(['..', filesep, 'fn_el_', new_el_type, '.m'], K, C, M, detJ, loc_nd, loc_df, no_dims);

%--------------------------------------------------------------------------
%TEST - SINGLE ELEMENT
no_nds = size(nds_in_nat_coords, 1);
test_nds = rand(no_nds, no_dims);
test_els = 1:no_nds;

fn_el_mats_test = str2func(['fn_el_', new_el_type]);
fn_el_mats_ref = str2func(['fn_el_', ref_el_type]);

test_D = rand(6);
test_D = test_D + test_D';
test_rho = 1234.5;

%New function
[test_el_K, test_el_C, test_el_M, test_loc_nd, test_loc_df] = fn_el_mats_test(test_nds, test_els, test_D, test_rho, [1,2,3]);

if exist(func2str(fn_el_mats_ref), 'file')
    %Existing function
    [el_K, el_C, el_M, loc_nd, loc_df] = fn_el_mats_ref(test_nds, test_els, test_D, test_rho, [1,2,3]);
    %Comparison
    fprintf('\nCOMPARISON WITH EXISTING\n')
    fprintf('Fractional RMS error for K: %e\n', fn_compare_matrices(test_el_K, el_K));
    fprintf('Fractional RMS error for M: %e\n', fn_compare_matrices(test_el_M, el_M));
end

%--------------------------------------------------------------------------
%TEST - MANY ELEMENTS SPEED COMPARISON

fprintf('\nSPEED COMPARISON WITH EXISTING\n')
n = 1000000;
test_nds = rand(n, no_dims);
test_els = randi(n, n, no_nds);
if exist(func2str(fn_el_mats_ref), 'file')
    tic;
    [el_K2, el_C2, el_M2, loc_nd2, loc_df2] = fn_el_mats_ref(test_nds, test_els, test_D, test_rho);
    fprintf('Existing form took %.2fs\n', double(toc));
end
tic;
[el_K2, el_C2, el_M2, loc_nd2, loc_df2] = fn_el_mats_test(test_nds, test_els, test_D, test_rho);
fprintf('New form took %.2fs\n', double(toc));

%--------------------------------------------------------------------------
%TESTS - AGAINST POGO (needs small model - Pogo cannot handle single
%element model)

fe_options.pogo_path = 'C:\Program Files\Pogo\windows\new version';
fe_options.pogo_matlab_path = 'C:\Program Files\Pogo\matlab';


%Say which element type will be used (currently there is only one choice for
%a solid material in 2D model, but in the future there may be more options to choose
%from, e.g. quadrilateral elements, second order elements)
el_typ_to_use_for_solid = new_el_type; 

model_size = 0.01;
bdry_pts = [
    0,          0 
    model_size, 0 
    model_size, model_size
    0,          model_size];

%Define the material
matls{1} = fn_matl_isotropic_solid_defined_by_velocities('steel', 5800, 3000, 8900);

%Work out element size and time step
el_size = 0.01;
time_step = fn_get_suitable_time_step(matls, el_size);
max_time = 10 * time_step; 

%Create the nodes and elements of the mesh
mod = fn_2d_structured_mesh_rectangular_els(bdry_pts, el_size);

%Associate each element with a material index and element type index
mod.el_mat_i(:) = 1;
mod.el_typ_i(:) = 1;

%Identify node closest to desired source location
steps{1}.load.frc_nds = 1;
steps{1}.load.frc_dfs = 1;

%Provide the time signal for the loading
steps{1}.load.time = 0: time_step:  max_time;
steps{1}.load.frcs = zeros(size(steps{1}.load.time));
steps{1}.load.frcs(1) = 1;

%Say where the displacement should be monitored
steps{1}.mon.nds = 1;
steps{1}.mon.dfs =  1;

solvers = {'BristolFE', 'pogo'};
fe_options.dof_to_use = [1,2];
fe_options.solver_precision = 'double';
for s = 1:numel(solvers)
    fe_options.solver = solvers{s};
    switch solvers{s}
        case 'BristolFE'
            el_types = {new_el_type};
        case 'pogo'
            el_types = {ref_el_type};
    end
    [~, tmp] = fn_FE_entry_point(mod, matls, el_types, steps, fe_options);
    switch solvers{s}
        case 'BristolFE'
            test_mod_K = tmp.K;
            test_mod_M = tmp.M;
        case 'pogo'
            pogo_mod_K = tmp.K;
            pogo_mod_M = tmp.M;
    end

end

fprintf('\nCOMPARISON WITH POGO\n')
fprintf('Fractional RMS error for model K compared to pogo: %e\n', full(fn_compare_matrices(test_mod_K, pogo_mod_K)));
fprintf('Fractional RMS error for model M compared to pogo: %e\n', full(fn_compare_matrices(test_mod_M, pogo_mod_M)));
