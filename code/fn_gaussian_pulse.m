function s = fn_gaussian_pulse(t, centre_freq, no_cycles, varargin)
%USAGE
%   s = fn_gaussian_pulse(t, centre_freq, no_cycles[, db_down_at_start, db_down])
%AUTHOR
%   Paul Wilcox (2025)
%SUMMARY
%   Returns a Gaussian pulse with specifed centre frequency and number of
%   cycles given a specified time axis.
%INPUTS
%   t - vector of time points
%   centre_freq - centre frequency of pulse
%   no_cycles - number of cycles in pulse (based on -db_down points)
%   [db_down_at_start - amplitude of envelope at t(1), default = 60]
%   [db_down - how length of pulse is defined, default = 40]
%OUTPUTS
%   s - the generated Gaussian pulse
%NOTES
%   Pulse will start from start of time-axis. First point of s will be 
%   first point where envelope of pulse > -20*log10(db_down_at_start). 
%   Output will have same size as t.
%--------------------------------------------------------------------------

if numel(varargin) >= 1 && ~isempty(varargin{1})
    db_down_at_start = abs(varargin{1});
else
    db_down_at_start = 60;
end
if numel(varargin) >= 2 && ~isempty(varargin{2})
    db_down = abs(varargin{2});
else
    db_down = 40;
end

sz = size(t);
t = t(:);
beta = 10 ^ (-db_down / 20);
T = no_cycles / (2 * centre_freq * sqrt(log(1/beta)));
gamma = 10 ^ (-db_down_at_start / 20);
ct = t(1) + sqrt(-T ^ 2 * log(gamma));

env = exp(-((t - ct)/T) .^ 2);
s = env .* sin(2*pi * centre_freq * (t - ct));
s = reshape(s, sz);
end