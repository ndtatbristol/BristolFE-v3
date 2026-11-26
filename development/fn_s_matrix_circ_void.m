function S = fn_s_matrix_circ_void(inc_angles_degs, k, a)
ka = k * a;
theta = inc_angles_degs * pi / 180;
nmax = 20;
S = zeros(size(inc_angles_degs));
for n = -nmax:nmax
    S = S + ...
        (n / ka * besselj(n, ka) - besselj(n + 1, ka)) ./ ...
        (n / ka * besselh(n, 1, ka) - besselh(n + 1, 1, ka)) * exp(1i * n * theta);
end
S = -S * sqrt(2*pi/k) / (pi * exp(1i * pi / 4));
end