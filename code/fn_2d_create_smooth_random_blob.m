function pts = fn_2d_create_smooth_random_blob(min_rad_frac, complexity, no_pts)
%USAGE
%   pts = fn_2d_create_smooth_random_blob(min_rad_frac, complexity, no_pts)
%AUTHOR
%   Paul Wilcox (2025)
%SUMMARY
%   Creates list of 2D points that describe perimeter of smooth blob
%   shape. 
%INPUTS
%   min_rad_frac - minimum radial size
%   complexity - number that determines how many orders of random harmonics
%   are used to generate shape, hence complexity = 0 just gives a circle of
%   radius 1.
%   no_pts - how many points in output
%OUTPUT
%   pts - no_pts x 2 matris of perimeter coordinates of blob. Centre is 
%   always at (0,0) and maximum radial distance of any perimeter point is 
%   1, so result should be scaled by desired max radius and coordinates of
%   desired centre added on, e.g.
%       actual_pts = pts * desired_max_radius + desired_centre
%NOTES
%   min_rad_frac = 0.5 and complexity = 3 give nice smooth blobby shapes
%--------------------------------------------------------------------------
a = linspace(0, 2 * pi, no_pts + 1)';
a = a(1:end - 1);
if complexity > 0
    n = [0: complexity];
    phi = rand(size(n)) * 2 * pi;
    r = sum(sin(a * n + phi), 2);
    r = (r - min(r)) / (max(r) - min(r)) * (1 - min_rad_frac) + min_rad_frac;
else
    r = 1;
end
pts = r .* [cos(a), sin(a)];
end
