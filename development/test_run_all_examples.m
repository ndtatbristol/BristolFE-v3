clear all;
close all;

%ABOUT THIS SCRIPT
%This runs all of the example models, with low number of elements per 
%wavelength and shows results. Intendeded for quick check that all other
%exampes will work.

%Uncomment one of the following model file names to determine which one
%will be used for comparison:
models_to_run = {@mod_2d_basic, @mod_3d_basic, @mod_2d_advanced, @mod_3d_advanced, @mod_2d_oblique};
solvers_to_use = {'BristolFE', 'pogo'};

% models_to_run = {@mod_2d_basic, @mod_3d_basic};

%Parameters for all models
fixed_params.els_per_wavelength = 3;
fixed_params.include_fluid_region = 0;
fixed_params.fe_options_field_output_every_n_frames = inf;
fixed_params.fe_options_solver_mode = 'pc';
fixed_params.fe_options_sort_nds = 1;
fixed_params.fe_options_pogo_path = 'C:\Program Files\pogo\bin'; %Ignored if solver is not pogo
fixed_params.fe_options_pogo_matlab_path = 'C:\Program Files\Pogo\matlab'; %Ignored if solver is not pogo

%--------------------------------------------------------------------------

%Add all Bristol FE functions to Matlab path
addpath(genpath([fileparts(mfilename('fullpath')), filesep, '..', filesep, 'code']))
addpath(genpath([fileparts(mfilename('fullpath')), filesep, '..', filesep, 'examples', filesep, 'models']))

%Define the model
for m = 1:numel(models_to_run)
    for s = 1:numel(solvers_to_use)
        params = fixed_params;
        params.fe_options_solver = solvers_to_use{s};
        model_to_run = models_to_run{m};
        str = ['Model: ', func2str(model_to_run), '; Solver: ' , params.fe_options_solver];
        fprintf('\n---------------------------------------------------------------------------\n')
        fprintf(str);
        fprintf('\n---------------------------------------------------------------------------\n')

        %Define the model
        [mod, matls, el_types, steps, fe_options, params] = model_to_run(params);
    
        %Execute the model
        res = fn_FE_entry_point(mod, matls, el_types, steps, fe_options);
        
        %Plot summed history output over monitoring nodes
        figure;
        plot(steps{1}.load.time, res{1}.dsps);
        % plot(steps{1}.load.time, sum(res{1}.dsps));
        xlabel('Time (s)')
        title(str, 'Interpreter', 'none')
    end
end

