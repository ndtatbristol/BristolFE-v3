%test new element formulation method
clear all
addpath('../')

%CPE3
nds_in_nat_coords = [
    0, 0; 
    1, 0; 
    0, 1];
sf_powers = [
    0, 0
    1, 0
    0, 1];
gauss_pts = [
    1/3, 1/3];
gauss_weights = 1/2; 


%CPE4
% nds_in_nat_coords = [
%     -1, -1 
%      1, -1 
%      1,  1 
%     -1,  1];
% sf_powers = [
%     0, 0 
%     1, 0 
%     0, 1 
%     1, 1];

no_dfs = 3;
solid_or_fluid = 'solid';
el_type = 'test';
no_dims = 2;


[K, loc_nd, loc_df] = fn_symbolic_K_matrix(nds_in_nat_coords, gauss_pts, gauss_weights, sf_powers, no_dfs, solid_or_fluid);

C = sym(zeros(size(K)));
M = sym(zeros(size(K)));
detJ = [];
%--------------------------------------------------------------------------
fn_create_element_matrix_file(['..', filesep, 'fn_el_', el_type, '.m'], K, C, M, detJ, loc_nd, loc_df, no_dims);

fn_el_mats_test = str2func(['fn_el_', el_type]);
D = rand(6);
D = D + D';
rho = 1234.5;

%Test function #1 - limited DOF, just one element
[el_K, el_C, el_M, loc_nd, loc_df] = fn_el_mats_test([0,0;1,0;0,1], [1,2,3], D, rho, [1,2]);
disp(squeeze(el_K));
[el_K, el_C, el_M, loc_nd, loc_df] = fn_el_CPE3([0,0;1,0;0,1], [1,2,3], D, rho, [1,2]);
disp(squeeze(el_K));
