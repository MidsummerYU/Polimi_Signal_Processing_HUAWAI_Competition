%% Solutions

close all; clear sound; clear; clc;
load('Hw2a.mat'); load('Hw2b.mat'); load('Hw2c.mat');
load('Hw2d.mat'); load('Hw2e.mat');




%% Listen to Original Sin_a
sound(Sin_a, fs);

%% Stop
clear sound;

%% Create, Plot & Play S_est_a
S_est_a(:,1) = WLMS(Sin_a(:,1), Sn_ref_a(:,1), 0.06, 60);
S_est_a(:,2) = WLMS(Sin_a(:,2), Sn_ref_a(:,2), 0.06, 60);

S_est_a(:,1) = hardcut_filter(S_est_a(:,1), fs, 3300, 5200);
S_est_a(:,2) = hardcut_filter(S_est_a(:,2), fs, 3300, 5200);

S_est_a(:,1) = S_est_a(:,1) / max(abs(S_est_a(:,1)) * 1.1);
S_est_a(:,2) = S_est_a(:,2) / max(abs(S_est_a(:,2)) * 1.1);

sound(S_est_a, fs);
figure;
spectrogram(S_est_a(:,1), hamming(2048), 1024, 2048, fs, 'yaxis');
colorbar;
title('Spectrogram of S_{est a}');

%% Stop
clear sound;

%% Store
to_store = S_est_a ./ max(abs(S_est_a(:))) * 0.99;
audiowrite('HW2_a.mp3', to_store, fs);




%% Listen to Original Sin_b
sound(Sin_b, fs);

%% Stop
clear sound;

%% Create, Plot & Play S_est_b
S_est_b(:,1) = SWLMS(Sin_b(:,1), Sn_ref_b(:,2), 0.03, 500, fs);
S_est_b(:,2) = SWLMS(Sin_b(:,2), Sn_ref_b(:,1), 0.03, 500, fs);

S_est_b(:,1) = S_est_b(:,1) / max(abs(S_est_b(:,1)) * 1.1);
S_est_b(:,2) = S_est_b(:,2) / max(abs(S_est_b(:,2)) * 1.1);

sound(S_est_b, fs);
figure;
spectrogram(S_est_b(:,1), hamming(2048), 1024, 2048, fs, 'yaxis');
colorbar;
title('Spectrogram of S_{est b}');

%% Stop
clear sound;

%% Store
to_store = S_est_b ./ max(abs(S_est_b(:))) * 0.99;
audiowrite('HW2_b.mp3', to_store, fs);




%% Listen to Original Sin_c
sound(Sin_c, fs);

%% Stop
clear sound;

%% Create, Plot & Play S_est_c
S_est_c(:,1) = SWLMS(Sin_c(:,1), Sn_ref_c(:,2), 0.03, 500, fs);
S_est_c(:,2) = SWLMS(Sin_c(:,2), Sn_ref_c(:,1), 0.03, 500, fs);

S_est_c(:,1) = S_est_c(:,1) / max(abs(S_est_c(:,1)) * 1.1);
S_est_c(:,2) = S_est_c(:,2) / max(abs(S_est_c(:,2)) * 1.1);

sound(S_est_c, fs);
figure;
spectrogram(S_est_c(:,1), hamming(2048), 1024, 2048, fs, 'yaxis');
colorbar;
title('Spectrogram of S_{est c}');

%% Stop
clear sound;

%% Store
to_store = S_est_c ./ max(abs(S_est_c(:))) * 0.99;
audiowrite('HW2_c.mp3', to_store, fs);




%% Listen to Original Sin_d
sound(Sin_d, fs);

%% Stop
clear sound;

%% Create, Plot & Play S_est_d
sin1 = create_stable_sin(1940, 816000, 11025, 0.02);
sin2 = create_stable_sin(1950, 408000, 11025, 0.02);
sin3 = create_variable_sin(1350, 500, 816000, 11025, 0.02);
sin4 = create_variable_sin(550, 500, 816000, 11025, 0.02);
sin5 = create_variable_sin(1950, 750, 816000, 11025, 0.02);
sin5 = sin5(1:408000);
Sn_ref_d = sin1+[sin2; sin5]+sin3+sin4;

S_est_d(:,1) = WLMS(Sin_d(:,1), Sn_ref_d, 0.036, 42);
S_est_d(:,2) = WLMS(Sin_d(:,2), Sn_ref_d, 0.049, 25);
S_est_d(:,1) = S_est_d(:,1) / (max(abs(S_est_d(:,1))) * 1.1);
S_est_d(:,2) = S_est_d(:,2) / (max(abs(S_est_d(:,2))) * 1.1);

sound(S_est_d, fs);
figure;
spectrogram(S_est_d(:,1), hamming(2048), 1024, 2048, fs, 'yaxis');
colorbar;
title('Spectrogram of S_{est d}');

%% Stop
clear sound;

%% Store
to_store = S_est_d ./ max(abs(S_est_d(:))) * 0.99;
audiowrite('HW2_d.mp3', to_store, fs);




%% Listen to Original Sin_e
sound(Sin_e, fs);

%% Stop
clear sound;

%% Create, Plot & Play S_est_e
Sn_ref_e = Sin_e - Sin_d;

S_est_e(:,1) = WLMS(Sin_e(:,1), Sn_ref_e(:,1), 0.008, 900);
S_est_e(:,2) = WLMS(Sin_e(:,2), Sn_ref_e(:,2), 0.008, 900);
S_est_e(:,1) = S_est_e(:,1) / (max(abs(S_est_e(:,1))) * 1.1);
S_est_e(:,2) = S_est_e(:,2) / (max(abs(S_est_e(:,2))) * 1.1);

figure;
spectrogram(S_est_e(:,1), hamming(2048), 1024, 2048, fs, 'yaxis');
colorbar;
title('Spectrogram of S_{est e}');
sound(S_est_e, fs);

%% Stop
clear sound;

%% Store
to_store = S_est_e ./ max(abs(S_est_e(:))) * 0.99;
audiowrite('HW2_e.mp3', to_store, fs);







%% Functions

% Creates a sine wave with constant frequency
function ref_signal = create_stable_sin(frequency, num_samples, fs, reference_amplitude)
    t = (0:num_samples-1)'/fs;
    ref_signal = sin(2*pi*frequency*t);
    scale_factor = rms(reference_amplitude)/rms(ref_signal);
    ref_signal = ref_signal * scale_factor;
end

% Creates a sine wave with frequency varying as a sine wave
function ref_signal = create_variable_sin(center_freq, freq_deviation, num_samples, fs, reference_amplitude)
    t = (0:num_samples-1)'/fs;
    mod_freq = fs/num_samples;  % Frequency for exactly one period
    inst_freq = center_freq + freq_deviation * sin(2*pi*mod_freq*t);
    phase = 2*pi * cumsum(inst_freq)/fs;
    ref_signal = sin(phase);
    scale_factor = rms(reference_amplitude)/rms(ref_signal);
    ref_signal = ref_signal * scale_factor;
end

% Computes Whitened LMS
function S_est = WLMS(x, ref, mu, filter_length)

    % Compute pre-whitening filter
    whitening_order = 5;
    [a_white, ~] = lpc(ref, whitening_order);
    a_white = a_white * 0.99;

    disp(length(a_white));
    
    % Pre-whiten reference signal
    ref_white = filter(a_white, 1, ref);
    
    % Normalize signals to prevent saturation
    ref_white = ref_white / max(abs(ref_white));
    x_normalized = x / max(abs(x));
    
    % LMS adaptation
    N = length(x);
    S_est = zeros(N, 1, 'like', x);
    w = zeros(filter_length, 1);
    ref_buffer = zeros(filter_length, 1);

    for n = 1:N
        ref_buffer = [ref_white(n); ref_buffer(1:end-1)];
        S_est(n) = x_normalized(n) - w' * ref_buffer;
        w = w + mu * S_est(n) * ref_buffer;
    end
    % Restore signal amplitude
    S_est = S_est * max(abs(x));
end

% Computes Subband Whitened LMS
function S_est = SWLMS(x, ref, step_size, filter_length, fs)
    N = length(x);
    
    % Define frequency bands
    bands = [
        0,    1000;  % Band 1
        1000, 2000;  % Band 2
        2000, 3000;  % Band 3
        3000, 4000;  % Band 4
        4000, fs/2;  % Band 5
    ];
    num_bands = size(bands, 1);
    
    % Initialize result
    S_est = zeros(N, 1);
    
    % Single loop for filtering and processing
    for i = 1:num_bands
        % Create appropriate filter based on band position
        if i == 1
            [b, a] = butter(6, bands(i,2)/(fs/2), 'low');
        elseif i == num_bands
            [b, a] = butter(6, bands(i,1)/(fs/2), 'high');
        else
            [b, a] = butter(6, [bands(i,1)/(fs/2), bands(i,2)/(fs/2)], 'bandpass');
        end
        
        % Filter signals and process immediately
        x_filtered = filtfilt(b, a, x);
        ref_filtered = filtfilt(b, a, ref);
        
        % Apply LMS and add to result directly
        S_est = S_est + WLMS(x_filtered, ref_filtered, step_size, round(filter_length));
    end
end

% Creates a simple filter to remove unwanted frequencies
function filtered_signal = hardcut_filter(signal, fs, f_low, f_high)
    % FFT of the signal
    N = length(signal);
    X = fft(signal);
    
    % Create filter mask
    mask = ones(N, 1);
    
    % Find indices corresponding to the cutoff frequencies
    idx_low = round(f_low*N/fs);
    idx_high = round(f_high*N/fs);
    
    % Set the mask to zero for the frequencies we want to cut
    mask(idx_low:idx_high) = 0;
    mask(N-idx_high:N-idx_low) = 0;  % Mirror for negative frequencies
    
    % Apply filter in frequency domain
    X_filtered = X .* mask;
    
    % Convert back to time domain
    filtered_signal = real(ifft(X_filtered));
end