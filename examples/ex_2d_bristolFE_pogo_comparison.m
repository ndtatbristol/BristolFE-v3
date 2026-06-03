clear all;
close all;

%ABOUT THIS SCRIPT
%This script runs the same basic model twice, once in BristolFE and once in Pogo
%(which must be installed with a valid license) and compares the results.

%Following will need to be set to where the Pogo executable and Pogo Matlab
%scripts are located respectively
pogo_path = 'C:\Program Files\Pogo\windows\new version';
pogo_matlab_path = 'C:\Program Files\Pogo\matlab';

params.element_shape = 'tri';

%--------------------------------------------------------------------------
%DEFINE THE MODEL

%Add all Bristol FE functions to Matlab path
addpath(genpath([fileparts(mfilename('fullpath')), filesep, '..', filesep, 'code']))

%Add models subfolder to path
addpath(['.', filesep, 'models']);

%Define the model
[mod, matls, el_types, steps, fe_options, params] = mod_2d_basic(params);

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
    plot(steps{1}.load.time, res{s}{1}.dsps, cols{s});
    hold on;
    xlabel('Time (s)')
end
legend(solvers)
subplot(2,1,2);
plot(steps{1}.load.time, 20*log10(abs(res{1}{1}.dsps - res{2}{1}.dsps) / max(abs(res{1}{1}.dsps))));
xlabel('Time (s)')
ylabel('Difference (db)')

