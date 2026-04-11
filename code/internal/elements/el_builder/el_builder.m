%General script for building and testing elements

clear all
rng(1);
cd(fileparts(mfilename('fullpath')));
addpath(genpath('..\..\..\..\code'));

simplify_expression = 0;

%CPE3 - OK
solid_or_fluid = 'solid';
new_el_type = 'CPE3_new';
ref_el_type = 'CPE3';
[nds_in_nat_coords, sf_powers, gauss_pts, gauss_weights, no_dims] = fn_el_parent_nds_and_shape_functions('triangular');

%CPE4 - OK
% solid_or_fluid = 'solid';
% new_el_type = 'CPE4_new4';
% ref_el_type = 'CPE4_new';
% [nds_in_nat_coords, sf_powers, gauss_pts, gauss_weights, no_dims] = fn_el_parent_nds_and_shape_functions('quadrilateral');

%AC2D3 - OK
% solid_or_fluid = 'fluid';
% new_el_type = 'AC2D3_new';
% ref_el_type = 'AC2D3';
% [nds_in_nat_coords, sf_powers, gauss_pts, gauss_weights, no_dims] = fn_el_parent_nds_and_shape_functions('triangular');

%AC2D4 - 
% solid_or_fluid = 'fluid';
% new_el_type = 'AC2D4_new';
% ref_el_type = '';
% [nds_in_nat_coords, sf_powers, gauss_pts, gauss_weights, no_dims] = fn_el_parent_nds_and_shape_functions('quadrilateral');


%C3D8
% solid_or_fluid = 'solid';
% new_el_type = 'C3D8_new3';
% ref_el_type = 'C3D8_new';
% [nds_in_nat_coords, sf_powers, gauss_pts, gauss_weights, no_dims] = fn_el_parent_nds_and_shape_functions('hexahedral');


%--------------------------------------------------------------------------
%DO THE SYMBOLIC CALCULATIONS AND CREATE THE ELEMENT FILE
new_el_fname = ['..', filesep, 'fn_el_', new_el_type, '.m'];

factorisation_level = 3;

%Symbolic calculation
sym_mats = fn_element_symbolic_matrices(nds_in_nat_coords, gauss_pts, gauss_weights, sf_powers, solid_or_fluid, factorisation_level);

%Write the element file
fn_create_element_matrix_file3(new_el_fname, sym_mats);

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
