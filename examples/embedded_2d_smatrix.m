clear
close all

%This is work in progress!!

model_to_run = @mod_2d_embedded_smatrix;

%Parameters for the model - if empty, default values for all parameters 
%will be used
params = [];

%However, any of the default parameters (see top of model file for complete 
%list in each case) can be overwritten here, e.g.
params.els_per_wavelength = 5;

%If you just want to see the model (without running it, set show_geom_only to 1
show_geom_only = 0;

%--------------------------------------------------------------------------
%DEFINE THE MODEL

%Add all Bristol FE functions to Matlab path
addpath(genpath([fileparts(mfilename('fullpath')), filesep, '..', filesep, 'code']))
addpath(genpath([fileparts(mfilename('fullpath')), filesep, '..', filesep, 'subdoms']))

%Add models subfolder to path
addpath(['.', filesep, 'models']);

%Define the model - could actually just return subdomain_mod as the only
%one since main is not required here
[mod, matls, el_types, steps, fe_options, params, subdomain_mod] = model_to_run(params);

%Need to modify above so that steps returns the displacements required at
%BLs 2 and 3 (and monitoring in BLs 1-4). These need to act as the recorded
%displacements from a main model so that usual subdomain code can be run

%Show the subdomain mesh and stop if requested
if show_geom_only
    figure;
    col = 'rgbc';
    display_options = [];
    for i = 1:4
        display_options.node_sets_to_plot(i).nd = find(subdomain_mod.bdry_lyrs == i);
        display_options.node_sets_to_plot(i).col = [col(i), '.'];
    end
    h_patch = fn_show_geometry(subdomain_mod, matls, el_types, display_options);
    return
end

main.doms{1}.mod = subdomain_mod;

%Get the matrices from the global modal
fe_options.return_matrices_only = 1;
main.res.mats = fn_FE_entry_point(mod, matls, el_types, [], fe_options);



main = fn_run_subdomain_model(main, fe_options);

%use create sub-domain function as usual and add scatterer

%rather than running main model to generate incident field on boundary
%nodes, use analytical solution for incident plane wave from each incident
%direction (may need some tricks here to deal with time shifted delta
%functions). Will need to have specify time-shift back from sub-domain
%centre so that incident field hits boundary at t = 0 or a bit later.

%run sub-domain model as usual for each incident direction and project into
%each scattered direction as usual (scatterer-free sub-domain should giev
%near-zero but it will be less accurate than pure FE version, especially at
%higher frequencies)

%FFT responses, remove in and out phase-shift due to propagation delay and that should be it.

