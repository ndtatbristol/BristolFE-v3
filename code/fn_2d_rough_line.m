function pts = fn_2d_rough_line(npts, length, angle, rms_roughness, corr_len)
%USAGE
%   pts = fn_2d_rough_line(length, angle, rms_roughness, corr_len)
%AUTHOR
%   Paul Wilcox (2026)
%SUMMARY
%   Returns rough line relative to mean straight line described by rms 
%   roughness and correlation length
%INPUTS
%   npts - number of points in output
%   length - length of mean line
%   angle - angle of mean line
%   rms_roughness - RMS roughness
%   corr_len - correlation length
%OUTPUT
%   pts - npts x 2 matrix of coordinates on the rough line with first
%   point always (0, 0)
%NOTES
%   Angles are in radians. 
%--------------------------------------------------------------------------
u = linspace(0, length, npts);
v = randn(size(u));
k = 2 * pi * [0:npts - 1] / (npts * (u(2) - u(1)));
V = fft(v) .* exp(-(k * corr_len / 2) .^ 2);
v = real(ifft(V));
v = v / sqrt(mean(v .^ 2)) * rms_roughness;

m = [cos(angle), -sin(angle); sin(angle), cos(angle)];
pts = (m * [u; v])';
pts = pts - pts(1, :);
end