function [max_damping, damping_power_law, max_stiffness_reduction] = fn_optimum_absorbing_bdry_properties(abs_bdry_thickness, matls, centre_freq)
%USAGE
%   [max_damping, damping_power_law, max_stiffness_reduction] = fn_optimum_absorbing_bdry_properties(abs_bdry_thickness, matls, centre_freq)
%AUTHOR
%   Paul Wilcox (2025)
%SUMMARY
%   Returns what appear to be good (not necessarily optimum!) values for
%   the parameters of absorbing boundary layers
%INPUTS
%   abs_bdry_thickness - thickness of the absorbing boundary layers,
%   typically using 1-2 times the longest wavelength possible is a good bet
%   matls - cell array of materials used in mod
%   centre_frequency - cemtre frequency of excitation
%OUTPUTS
%   max_damping, damping_power_law, max_stiffness_reduction - the
%   parameters to be used in the fe_options structure when running the
%   model
%NOTES
%   Absorbing boundaries are never perfect and determining the optimum 
%   absorbing properties is an unsolved problem. This function
%   returns some parameters that seem to work OK for most cases based on
%   emprical tests, but it has not been exhaustively tested for all 
%   possible cases.
%--------------------------------------------------------------------------
max_lambda = fn_estimate_max_min_vels(matls) / centre_freq;

abs_bdry_thickness_lambda = abs_bdry_thickness / max_lambda;

p =  [6.7368  -13.4737   16.4211];


max_damping = polyval(p, abs_bdry_thickness_lambda) * centre_freq;
damping_power_law = 3;
max_stiffness_reduction = 0.01;

end