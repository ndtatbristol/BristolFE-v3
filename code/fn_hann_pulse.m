function s = fn_hann_pulse(t, centre_freq, no_cycles)
%USAGE
%   s = fn_hann_pulse(t, centre_freq, no_cycles)
%AUTHOR
%   Paul Wilcox (2025)
%SUMMARY
%   Returns a Hann-windowed pulse with specifed centre frequency and number of
%   cycles given a specified time axis.
%INPUTS
%   t - vector of time points
%   centre_freq - centre frequency of pulse
%   no_cycles - number of cycles in pulse (based on -db_down points)
%OUTPUTS
%   s - the generated Hann-windowed pulse
%NOTES
%   Pulse will start from start of time-axis. First point will be zero. 
%   Output will have same size as t.
%--------------------------------------------------------------------------

sz = size(t);
t = t(:);
T = no_cycles / centre_freq;
ct = t(1) + T / 2;

env = 0.5 * (1 + cos((t - ct) / T * 2 * pi)) .* (abs(t - ct) <= (T / 2));
s = env .* sin(2 * pi * centre_freq * (t - ct));
s = reshape(s, sz);
end