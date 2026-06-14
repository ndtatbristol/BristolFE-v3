function h = fn_convolve(f, g, dim, options)
default_options.convole_use_gpu_if_available = 1;
default_options.convolve_chunk_size = 1024;

options = fn_set_default_fields(options, default_options);
%Convolves f by g along specified dim

chunk_size  = options.convolve_chunk_size;

use_gpu = 0;
if options.convole_use_gpu_if_available
    gpu_present = fn_test_if_gpu_present_and_working;
    if gpu_present
	    use_gpu = 1;
        reset(gpuDevice);
    end
end

fn_console_output(sprintf('Convolving %i signals with %i pts (GPU = %i) ... ', prod(size(f)) / size(f,dim), size(f,dim), use_gpu));
t1 = clock;

fft_pts = size(f, dim) * 2;

% if use_gpu
%     f = gpuArray(f);
%     g = gpuArray(g);
% end
% 
% F = fft(f, fft_pts, dim);
% G = fft(g, fft_pts, dim);
% H = F .* G;
% 
% h = ifft(H, fft_pts, dim);
% h = h(1: size(f,1), 1: size(f,2));
% 
% if use_gpu
%     h = gather(h);
%     reset(gpuDevice);
% end
%-------------------------------
    % ---- sizes ----
    [nf1, nf2] = size(f);
    [ng1, ng2] = size(g);

    if dim == 1
        % FFT down rows, broadcasting across columns
        if ~(nf2 == ng2 || nf2 == 1 || ng2 == 1)
            error('For dim == 1, size(f,2) and size(g,2) must match or one must be 1.');
        end

        n_keep = nf1;              % crop back to size(f,1)
        n_out_cols = max(nf2, ng2);

        h = complex(zeros(n_keep, n_out_cols, 'like', f));

        for c0 = 1:chunk_size:n_out_cols
            c1 = min(c0 + chunk_size - 1, n_out_cols);
            cols = c0:c1;

            % Map output columns back to f and g columns
            if nf2 == 1
                f_cols = ones(size(cols));
            else
                f_cols = cols;
            end

            if ng2 == 1
                g_cols = ones(size(cols));
            else
                g_cols = cols;
            end

            f_chunk = f(:, f_cols);
            g_chunk = g(:, g_cols);

            if use_gpu
                f_chunk = gpuArray(f_chunk);
                g_chunk = gpuArray(g_chunk);
            end

            F = fft(f_chunk, fft_pts, 1);
            G = fft(g_chunk, fft_pts, 1);
            H = F .* G;

            h_chunk = ifft(H, fft_pts, 1);
            h_chunk = h_chunk(1:n_keep, :);

            if use_gpu
                h_chunk = gather(h_chunk);
            end

            h(:, cols) = h_chunk;
        end

    else
        % dim == 2
        % FFT across columns, broadcasting across rows
        if ~(nf1 == ng1 || nf1 == 1 || ng1 == 1)
            error('For dim == 2, size(f,1) and size(g,1) must match or one must be 1.');
        end

        n_keep = nf2;              % crop back to size(f,2)
        n_out_rows = max(nf1, ng1);

        h = complex(zeros(n_out_rows, n_keep, 'like', f));

        for r0 = 1:chunk_size:n_out_rows
            r1 = min(r0 + chunk_size - 1, n_out_rows);
            rows = r0:r1;

            % Map output rows back to f and g rows
            if nf1 == 1
                f_rows = ones(size(rows));
            else
                f_rows = rows;
            end

            if ng1 == 1
                g_rows = ones(size(rows));
            else
                g_rows = rows;
            end

            f_chunk = f(f_rows, :);
            g_chunk = g(g_rows, :);

            if use_gpu
                f_chunk = gpuArray(f_chunk);
                g_chunk = gpuArray(g_chunk);
            end

            F = fft(f_chunk, fft_pts, 2);
            G = fft(g_chunk, fft_pts, 2);
            H = F .* G;

            h_chunk = ifft(H, fft_pts, 2);
            h_chunk = h_chunk(:, 1:n_keep);

            if use_gpu
                h_chunk = gather(h_chunk);
            end

            h(rows, :) = h_chunk;
        end
    end
%-------------------------------

fn_console_output(sprintf('completed in %.2f secs\n', etime(clock, t1)), [], 0);
end