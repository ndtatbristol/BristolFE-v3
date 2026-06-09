clear all
close all;

%ABOUT THIS SCRIPT
%This script runs the same model twice, once in BristolFE and once in Pogo
%(which must be installed with a valid license) and compares the results.

%Uncomment one of the following model file names to determine which one
%will be used for comparison:
model_to_run = @mod_2d_basic;
model_to_run = @mod_3d_basic;
% model_to_run = @mod_2d_advanced;

%Following will need to be set to where the Pogo executable and Pogo Matlab
%scripts are located respectively
pogo_path = 'C:\Program Files\Pogo\windows\new version';
pogo_matlab_path = 'C:\Program Files\Pogo\matlab';

%Parameters for the model - if empty, default values for all parameters 
%will be used
params = [];

%However, any of the default parameters (see top of model file for complete 
%list in each case) can be overwritten here, e.g.
params.els_per_wavelength = 10;

%If you just want to see the model (without running it, set 
%show_geom_only to 1
show_geom_only = 0;

%--------------------------------------------------------------------------
%DEFINE THE MODEL

%Add all Bristol FE functions to Matlab path
addpath(genpath([fileparts(mfilename('fullpath')), filesep, '..', filesep, 'code']))

%Add models subfolder to path
addpath(['.', filesep, 'models']);

%Define the model
params.include_fluid_region = 0; %Pogo does not support fluids so turn this off for any models that might include fluid regions
[mod, matls, el_types, steps, fe_options, params] = model_to_run(params);

%Show the mesh and stop if requested
if show_geom_only 
    figure;
    display_options.transparency = 0.5;
    display_options.draw_elements = 0;
    display_options.node_sets_to_plot(1).nd = steps{1}.load.frc_nds;
    display_options.node_sets_to_plot(1).col = 'r.';
    h_patch = fn_show_geometry(mod, matls, el_types, display_options);
    return
end

%--------------------------------------------------------------------------
%RUN THE MODEL FOR EACH SOLVER

solvers = {'BristolFE', 'pogo'};
fe_options.field_output_every_n_frames = inf; %Pogo field output currently not supported
fe_options.pogo_path = pogo_path;
fe_options.pogo_matlab_path = pogo_matlab_path;
for s = 1:numel(solvers)
    fe_options.solver = solvers{s};
    res{s} = fn_FE_entry_point(mod, matls, el_types, steps, fe_options);
end

%--------------------------------------------------------------------------
%SHOW THE RESULTS

%Plot history output at monitoring node for both solvers
figure;
subplot(2,1,1);
cols = {'r', 'k--'};
for s = 1:numel(solvers)
    ascan{s} = sum(res{s}{1}.dsps, 1);
    plot(steps{1}.load.time, ascan{s}, cols{s});
    hold on;
    xlabel('Time (s)')
end
legend(solvers)
subplot(2,1,2);
plot(steps{1}.load.time, 20*log10(abs(ascan{1} - ascan{2} ) / max(abs(ascan{1}))));
xlabel('Time (s)')
ylabel('Difference (db)')


