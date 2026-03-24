function rel_err = fn_compare_matrices(M1, M2)
rms_amp = sqrt((mean(abs(M1) .^ 2, 'all') + mean(abs(M2) .^ 2, 'all')) / 2);
rms_diff = sqrt(mean(abs(M1 - M2) .^ 2, 'all'));
rel_err = rms_diff / rms_amp;
end