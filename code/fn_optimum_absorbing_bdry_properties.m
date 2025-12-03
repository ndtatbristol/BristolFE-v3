function [max_damping, damping_power_law, max_stiffness_reduction] = fn_optimum_absorbing_bdry_properties(abs_bdry_thickness, matls, centre_freq)

max_lambda = fn_estimate_max_min_vels(matls) / centre_freq;

abs_bdry_thickness_lambda = abs_bdry_thickness / max_lambda;

p =  [6.7368  -13.4737   16.4211];


max_damping = polyval(p, abs_bdry_thickness_lambda) * centre_freq;
damping_power_law = 3;
max_stiffness_reduction = 0.01;

end