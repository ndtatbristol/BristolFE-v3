function [disp_output, force_output] = fn_freq_domain_solver_v1(...
    K, C, M, freq, ...
    forcing_indices, forcing_functions, ...
    disp_indices, disp_functions, varargin)

%SUMMARY
%   Solves dynamic problems in frequency domain. Output is always full
%   field as this is solved for anyway
%INPUTS
%   K - m x m global stiffness matrix
%   C - m x m global damping matrix
%   M - m x m global mass matrix
%   freq = n-element vector of frequencies
%   forcing_indices - p-element vector of global matrix indices at which 
%   forcing_functions will be applied (i.e. force input)
%   forcing_functions - p x s x n matrix of forces to apply
%   disp_indices - r-element vector of global matrix indices at which 
%   disp_functions will be applied (i.e. displacement input)
%   disp_functions - r x s x n matrix of displacments to apply
%OUTPUTS
%   disp_output - m x s x n matrix of displacements at all nodes
%   force_output - m x s x n matrix of external forces at all nodes (empty 
%   if no displacement input is used)

%--------------------------------------------------------------------------
if numel(varargin) < 1
	use_gpu_if_present = 1;
else
	use_gpu_if_present = varargin{1};
end
if numel(varargin) < 2
	solver_precision = 'double';
else
	solver_precision = varargin{2};
end

if isMATLABReleaseOlderThan('R2025a') && strcmpi(solver_precision, 'single')
    fprintf(['WARNING: Matlab before R2025a does not support sparse, single ' ...
        'precision matrices, defaulting to double precision instead\n']);
    solver_precision = 'double';
end

%--------------------------------------------------------------------------

gpu_present = fn_test_if_gpu_present_and_working;
if use_gpu_if_present && gpu_present
	use_gpu = 1;
else
	use_gpu = 0;
end

disp_control = 0;
load_control = 0;
%Error checks
if size(K, 1) ~= size(K, 2)
    error('K must be square matrix');
end
if size(C,1) ~= size(C, 2)
    error('C must be square matrix');
end
if size(M,1) ~= size(M, 2)
    error('M must be square matrix');
end
if ~isempty(forcing_indices)
    load_control = 1;
    if size(forcing_functions, 3) ~= numel(freq)
        error('size(forcing_functions, 3) does not equal number of frequencies');
    end
end
if ~isempty(disp_indices) 
    disp_control = 1;
    if size(disp_functions, 3) ~= numel(freq)
        error('size(disp_functions, 3) does not equal number of frequencies');
    end
end
if load_control && disp_control
    if size(forcing_functions, 2) ~= size(disp_functions, 2)
        error('Must be same number of forcing and displacement controlled load cases if both used');
    end
end

if ~isempty(forcing_indices)
    n_loads = size(forcing_functions, 2);
end
if ~isempty(disp_indices)
    n_loads = size(disp_functions, 2);
end
n_freqs = numel(freq);
n_dofs = size(K, 1);

fn_console_output(sprintf(['Frequency domain solver v1 (GPU = %i, load cases = %d, frequencies = %d, DoFs = %.3fM, ', solver_precision, ') ...'], use_gpu, n_loads, n_freqs, n_dofs / 1e6));

t1 = clock;
disp_output = zeros(n_dofs, n_loads, n_freqs);
frc = spalloc(n_dofs, n_loads, nnz(forcing_functions(:, :, 1)));
for fi = 1:n_freqs
    omega = freq(fi) * 2 * pi;
    X = K - 1i * omega * C - omega ^ 2 * M;
    frc(forcing_indices, :) = forcing_functions(:, :, fi);
    disp_output(:, :, fi) = X \ frc;
end

force_output = [];

fn_console_output(sprintf(' completed in %.2f secs\n', etime(clock, t1)), [], 0);

force_output = gather(force_output);
disp_output = gather(disp_output);
end