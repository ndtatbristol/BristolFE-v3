%test new element formulation method
clear all

%CPE3
nds_in_nat_coords = [
    0, 0; 
    1, 0; 
    0, 1];
sf_powers = [
    0, 0
    1, 0
    0, 1];

%CPE4
nds_in_nat_coords = [
    -1, -1 
     1, -1 
     1,  1 
    -1,  1];
sf_powers = [
    0, 0 
    1, 0 
    0, 1 
    1, 1];


shape_functions = fn_symbolic_shape_functions(nds_in_nat_coords, sf_powers);
N = fn_symbolic_shape_function_matrix(shape_functions, size(nds_in_nat_coords, 2));