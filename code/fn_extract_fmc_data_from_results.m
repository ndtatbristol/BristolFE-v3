function [time_data, tx, rx] = fn_extract_fmc_data_from_results(res)
%USAGE
%   [time_data, tx, rx] = fn_extract_fmc_data_from_results(res)
%AUTHOR
%   Paul Wilcox (2025)
%SUMMARY
%   Extracts results from FE model with multiple steps into single matrix 
%   of FMC data. 
%INPUTS
%   res - cell array output from multi-step FE model.
%OUTPUTS
%   time_data - 2D matrix of FMC data (consistent with Bristol UNDT group's
%   exp_data.time_data). Each column is an A-scan and transmission/reception
%   elements associated with each column are given in row vectors tx and rx
%   tx, tx - row vectors identifying transmitter and receiver assoicated
%   with each column in time_data
%--------------------------------------------------------------------------

%First pass just to work out how much data there is
m = 0;
n = 0;
for s = 1:numel(res)
    n = n + size(res{s}.dsps, 1);
    m = max(m, size(res{s}.dsps, 2));
end
time_data = zeros(m, n);
tx = zeros(1, n);
rx = zeros(1, n);

%second pass to populate
i1 = 1;
for s = 1:numel(res)
    i2 = i1 + size(res{s}.dsps, 1) - 1;
    tx(i1:i2) = s;
    rx(i1:i2) = (i1:i2) - i1 + 1;
    time_data(1:size(res{s}.dsps, 2), i1:i2) = res{s}.dsps.';
    i1 = i1 + size(res{s}.dsps, 1);
end

end


