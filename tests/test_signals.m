clear
close all

dt = 1/100e6;
t = [1:1000]' * dt;
centre_freq = 5e6;
no_cycles = 5;
db_down_at_start = 60;
db_down = 60;

gp = fn_gaussian_pulse(t, centre_freq, no_cycles, db_down_at_start, db_down);
hp = fn_hann_pulse(t, centre_freq, no_cycles);

n = numel(t);
df = 1 / (n * dt);
f = (0: n - 1)' * df;

figure;
subplot(1,2,1);
plot(t, [gp, hp]);
xlim([0, no_cycles / centre_freq * 2])
legend('Gaussian', 'Hann');

subplot(1,2,2);
plot(f, abs(fft([gp, hp], n)));
xlim([0, 3 * centre_freq]);
legend('Gaussian', 'Hann');
