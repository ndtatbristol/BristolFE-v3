%General script for building and testing elements

clear all
clc;

rng(1);
cd(fileparts(mfilename('fullpath')));
addpath(genpath('..\..\..\..\code'));

el_types = {'CPE3', 'AC2D3', 'CPE4', 'AC2D4'};
el_types = {'C3D8'};
ref_el_suffix = '_f3';
new_el_suffix = '_f3';
factorisation_level = 3;
no_trials = 100000;

build_element_functions = 1;
test_element_functions = 1;

if build_element_functions
    for el = 1:numel(el_types)
        el_type = el_types{el};
        
        %Get the details for the specified element type
        [nds_in_nat_coords, sf_powers, gauss_pts, gauss_weights, no_dims, solid_or_fluid] = fn_el_details_for_builder(el_type);

        new_el_fn_name = ['fn_el_', el_type, new_el_suffix];
        new_el_fname = ['..', filesep, new_el_fn_name, '.m'];

        %Symbolic calculation
        fprintf(['Generating symbolic matrices for ', el_type, ':\n']);
        sym_mats = fn_element_symbolic_matrices(nds_in_nat_coords, gauss_pts, gauss_weights, sf_powers, solid_or_fluid, factorisation_level);

        %Write the element file
        fprintf(['Writing element file for ', new_el_fn_name, '\n']);
        fn_write_element_matrix_file(new_el_fname, sym_mats);

        %Test it - this is purely to check for errors
        fprintf(['Performing numerical test of ', el_type, '\n']);
        [el_K, el_M, el_C, time_taken] = fn_test_element_numeric(new_el_fn_name, nds_in_nat_coords, 10);
    end
end

if test_element_functions
    for el = 1:numel(el_types)
        el_type = el_types{el};
        [nds_in_nat_coords, sf_powers, gauss_pts, gauss_weights, no_dims, solid_or_fluid] = fn_el_details_for_builder(el_type);        new_el_fn_name = ['fn_el_', el_type, new_el_suffix];
        ref_el_fn_name = ['fn_el_', el_type, ref_el_suffix];
        fprintf(['Testing ', new_el_fn_name, ' for %d elements:\n'], no_trials);
        [el_K1, el_M1, el_C1, time_taken1] = fn_test_element_numeric(new_el_fn_name, nds_in_nat_coords, no_trials);
        fprintf('\tTime taken: %.3f\n', time_taken1);
        if exist(ref_el_fn_name, 'file')
            fprintf(['Testing ', ref_el_fn_name, ' for %d elements:\n'], no_trials);
            [el_K2, el_M2, el_C2, time_taken2] = fn_test_element_numeric(ref_el_fn_name, nds_in_nat_coords, no_trials);
            fprintf('\tTime taken: %.3f\n', time_taken2);
            fprintf(['\nComparison between ', new_el_fn_name,' and ', ref_el_fn_name,':\n']);
            fprintf('  Fractional RMS difference for K: %e\n', fn_compare_matrices(el_K1, el_K2));
            fprintf('  Fractional RMS difference for M: %e\n', fn_compare_matrices(el_M1, el_M2));
            fprintf('  Fractional RMS difference for C: %e\n', fn_compare_matrices(el_C1, el_C2));
        end
    end
end
return
%--------------------------------------------------------------------------
%DO THE SYMBOLIC CALCULATIONS AND CREATE THE ELEMENT FILE



%Tests
test_D = rand(size(sym_mats.D));
test_D = test_D + test_D';
test_rho = rand(1);
test_nds = nds_in_nat_coords + randn(size(nds_in_nat_coords)) * 0.1;

%First test - just of the symbolic result (not the element file)
[K_test, M_test] = fn_test_symbolic_matrices(sym_mats, test_nds, test_D, test_rho);
[K_ref, M_ref] = fn_el_test(ref_el_type, '', solid_or_fluid, test_nds, test_D, test_rho, 1);
fprintf(['\nCOMPARISON OF OUTPUTS FROM NEW SYMBOLIC MATRICES (FACTORISATION %i) AND ', ref_el_type,':\n'], factorisation_level);
fprintf('  Fractional RMS error for K: %e\n', fn_compare_matrices(K_test, K_ref));
fprintf('  Fractional RMS error for M: %e\n', fn_compare_matrices(M_test, M_ref));

%Second test - test of new element file for one element
[K_test, M_test, K_ref, M_ref] = fn_el_test(new_el_type, ref_el_type, solid_or_fluid, test_nds, test_D, test_rho, 1);
fprintf(['\nCOMPARISON OF OUTPUTS FROM NEW IMPLEMENTATIONS (FACTORISATION %i) AND ', ref_el_type,':\n'], factorisation_level);
fprintf('  Fractional RMS error for K: %e\n', fn_compare_matrices(K_test, K_ref));
fprintf('  Fractional RMS error for M: %e\n', fn_compare_matrices(M_test, M_ref));

%Third test - speed test for lost of elements
fn_el_test(new_el_type, ref_el_type, solid_or_fluid, test_nds, test_D, test_rho, 10000);

return

%--------------------------------------------------------------------------
%This bit needs to be done better with proper models chosen according to
%element type

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
switch solid_or_fluid
    case 'solid'
        matls{1} = fn_matl_isotropic_solid_defined_by_velocities('steel', 5800, 3000, 8900);
        active_dof = 1;
        solvers = {'BristolFE', 'pogo'};
        fe_options.dof_to_use = 1:no_dims;
    case 'fluid'
        matls{1} = fn_matl_fluid_defined_by_velocity('water', 1480, 1000);
        active_dof = 4;
        solvers = {'BristolFE'};
        fe_options.dof_to_use = 4;
end

%Work out element size and time step
el_size = 0.01;
time_step = fn_get_suitable_time_step(matls, el_size);
max_time = 10 * time_step;

%Create the nodes and elements of the mesh
mod = fn_2d_structured_mesh_triangular_els(bdry_pts, el_size);

%Associate each element with a material index and element type index
mod.el_mat_i(:) = 1;
mod.el_typ_i(:) = 1;

%Identify node closest to desired source location
steps{1}.load.frc_nds = 1;
steps{1}.load.frc_dfs = active_dof;

%Provide the time signal for the loading
steps{1}.load.time = 0: time_step:  max_time;
steps{1}.load.frcs = zeros(size(steps{1}.load.time));
steps{1}.load.frcs(1) = 1;

steps{1}.mon.nds = 1;
steps{1}.mon.dfs = active_dof;

fe_options.solver_precision = 'double';
for s = 1:numel(solvers)
    fe_options.solver = solvers{s};
    switch solvers{s}
        case 'BristolFE'
            el_types = {new_el_type};
        case 'pogo'
            el_types = {ref_el_type};
    end
    [~, mats{s}] = fn_FE_entry_point(mod, matls, el_types, steps, fe_options);
end

if numel(solvers) > 1
    fprintf('\nCOMPARISON WITH POGO\n')
    fprintf('Fractional RMS error for model K compared to pogo: %e\n', full(fn_compare_matrices(mats{1}.K, mats{2}.K)));
    fprintf('Fractional RMS error for model M compared to pogo: %e\n', full(fn_compare_matrices(mats{1}.M, mats{2}.M)));
end
