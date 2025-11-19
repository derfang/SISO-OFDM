%% Project 1: Adaptive OFDM System (Final Consolidated)
% Objective: End-to-End Simulation of ZF vs MMSE with Adaptive Switching
% Features: Parallel Processing (BER) + Detailed Constellation Analysis
clear; clc; close all;

%% --- SECTION 1: SYSTEM PARAMETERS ---
N_fft     = 64;        % Subcarriers
N_cp      = 16;        % Cyclic Prefix
Mod_Order = 16;        % 16-QAM
k_bits    = log2(Mod_Order);

% Simulation Configuration
N_sym_BER = 5000;      % High symbol count for smooth BER curves (uses Parallel)
N_sym_Vis = 1000;      % Lower symbol count for Visualization (Single thread)
SNR_Range = 0:2:30;    % SNR range for Waterfall
SNR_Vis   = 20;        % Specific SNR for Constellation plots
Threshold = 15;        % Adaptive Switch Point (dB)

% Check for Parallel Computing Toolbox
if isempty(gcp('nocreate'))
    try
        parpool; % Start workers
    catch
        fprintf('Parallel Toolbox not found. Using single core.\n');
    end
end

fprintf('--- Starting Final Simulation ---\n');

%% --- SECTION 2: BER WATERFALL (Parallelized) ---
fprintf('1. Calculating BER Waterfall (Parallel)... \n');

BER_ZF       = zeros(length(SNR_Range), 1);
BER_MMSE     = zeros(length(SNR_Range), 1);
BER_Adaptive = zeros(length(SNR_Range), 1);
snr_lin_vec  = 10.^(SNR_Range./10);

parfor i = 1:length(SNR_Range)
    current_snr = SNR_Range(i);
    snr_lin     = snr_lin_vec(i);
    
    % Local counters for this worker
    local_err_ZF = 0; local_err_MMSE = 0; local_err_Adapt = 0;
    
    % Generate Data (Local to this thread)
    tx_bits_loc = randi([0 1], N_fft * N_sym_BER * k_bits, 1);
    tx_mod_loc  = qammod(tx_bits_loc, Mod_Order, 'InputType', 'bit', 'UnitAveragePower', true);
    X_tx_loc    = reshape(tx_mod_loc, N_fft, N_sym_BER);
    
    % Loop over symbols
    for sym_idx = 1:N_sym_BER
        % Tx
        X_freq = X_tx_loc(:, sym_idx);
        x_time = ifft(X_freq, N_fft);
        x_cp   = [x_time(end-N_cp+1:end); x_time];
        
        % Noise Calcs
        sig_pwr   = mean(abs(x_cp).^2);
        noise_pwr = sig_pwr / snr_lin;
        
        % Channel: Random Rayleigh Fading
        h = (randn(4,1) + 1j*randn(4,1)) / sqrt(2);
        h = h / norm(h);
        
        % Rx
        y_chan = conv(x_cp, h);
        y_chan = y_chan(1:length(x_cp)); % Truncate convolution tail
        
        noise = sqrt(noise_pwr/2) * (randn(size(y_chan)) + 1j*randn(size(y_chan)));
        y_rx  = y_chan + noise;
        
        % Equalizer Prep
        y_no_cp = y_rx(N_cp+1:end);
        Y_Rx    = fft(y_no_cp, N_fft);
        H_freq  = fft(h, N_fft);
        
        % Equalization
        % ZF
        X_hat_ZF = Y_Rx ./ H_freq;
        
        % MMSE
        N0_eff = noise_pwr * N_fft; 
        X_hat_MMSE = (conj(H_freq) .* Y_Rx) ./ (abs(H_freq).^2 + N0_eff);
        
        % Adaptive
        if current_snr < Threshold
            X_hat_Adapt = X_hat_MMSE;
        else
            X_hat_Adapt = X_hat_ZF;
        end
        
        % Demod & Count
        bits_ref = qamdemod(X_freq, Mod_Order, 'OutputType', 'bit', 'UnitAveragePower', true);
        
        bits_ZF    = qamdemod(X_hat_ZF, Mod_Order, 'OutputType', 'bit', 'UnitAveragePower', true);
        bits_MMSE  = qamdemod(X_hat_MMSE, Mod_Order, 'OutputType', 'bit', 'UnitAveragePower', true);
        bits_Adapt = qamdemod(X_hat_Adapt, Mod_Order, 'OutputType', 'bit', 'UnitAveragePower', true);

        local_err_ZF    = local_err_ZF    + biterr(bits_ref, bits_ZF);
        local_err_MMSE  = local_err_MMSE  + biterr(bits_ref, bits_MMSE);
        local_err_Adapt = local_err_Adapt + biterr(bits_ref, bits_Adapt);
    end
    
    % Store Results
    BER_ZF(i)       = local_err_ZF / (N_sym_BER * N_fft * k_bits);
    BER_MMSE(i)     = local_err_MMSE / (N_sym_BER * N_fft * k_bits);
    BER_Adaptive(i) = local_err_Adapt / (N_sym_BER * N_fft * k_bits);
    
    fprintf('Completed %d dB\n', current_snr);
end

%% --- SECTION 3: CONSTELLATION ANALYSIS (Single Thread) ---
fprintf('2. Generating Constellation Visuals (Fixed Channel)...\n');

% Fixed Channel for consistent "Donut" visual
h_vis = [1; 0.8; 0.3]; 
h_vis = h_vis / norm(h_vis);

% Data Setup
tx_bits_vis = randi([0 1], N_fft * N_sym_Vis * k_bits, 1);
tx_mod_vis  = qammod(tx_bits_vis, Mod_Order, 'InputType', 'bit', 'UnitAveragePower', true);
X_tx_vis    = reshape(tx_mod_vis, N_fft, N_sym_Vis);

% Tx
x_time_vis = ifft(X_tx_vis, N_fft);
x_cp_vis   = [x_time_vis(end-N_cp+1:end, :); x_time_vis];
x_vec_vis  = x_cp_vis(:);

% Noise Setup
sig_pwr_vis   = mean(abs(x_vec_vis).^2);
snr_lin_vis   = 10^(SNR_Vis/10);
noise_pwr_vis = sig_pwr_vis / snr_lin_vis;

% Rx (Correct Dimensions Fix)
y_chan_vis = conv(x_vec_vis, h_vis);
y_chan_vis = y_chan_vis(1:length(x_vec_vis)); % Truncate BEFORE adding noise

noise_vis  = sqrt(noise_pwr_vis/2) * (randn(size(y_chan_vis)) + 1j*randn(size(y_chan_vis)));
y_rx_vis   = y_chan_vis + noise_vis;

% Process
y_mat_vis = reshape(y_rx_vis, N_fft+N_cp, N_sym_Vis);
Y_Rx_Vis  = fft(y_mat_vis(N_cp+1:end, :), N_fft);

% Equalize
H_freq_vis = fft(h_vis, N_fft);
H_mat_vis  = repmat(H_freq_vis, 1, N_sym_Vis);
N0_vis     = noise_pwr_vis * N_fft;

Rx_ZF_Vis   = Y_Rx_Vis ./ H_mat_vis;
Rx_MMSE_Vis = (conj(H_mat_vis) .* Y_Rx_Vis) ./ (abs(H_mat_vis).^2 + N0_vis);

%% --- SECTION 4: FINAL PLOTTING ---

% Figure 1: The Science (BER Waterfall)
figure('Name', 'BER Performance', 'Color', 'w', 'Position', [100 100 800 600]);
semilogy(SNR_Range, BER_ZF, '-ro', 'LineWidth', 2, 'MarkerSize', 6);
hold on;
semilogy(SNR_Range, BER_MMSE, '-b^', 'LineWidth', 2, 'MarkerSize', 6);
semilogy(SNR_Range, BER_Adaptive, '--k', 'LineWidth', 2);
xline(Threshold, 'g--', 'Switch Threshold');
grid on;
legend('Zero Forcing', 'MMSE', 'Adaptive', 'Location', 'southwest');
xlabel('SNR (dB)'); ylabel('Bit Error Rate (BER)');
title('OFDM Performance: ZF vs MMSE');
axis([min(SNR_Range) max(SNR_Range) 1e-5 1]);

% Figure 2: The Story (Constellation Transformation)
figure('Name', 'Constellation Analysis', 'Color', 'w', 'Position', [150 150 1200 450]);
ideal_pts = qammod((0:15)', 16, 'UnitAveragePower', true);

% Plot 1: Raw
subplot(1, 3, 1);
plot(Y_Rx_Vis(:), '.', 'Color', [0.6 0.6 0.6], 'MarkerSize', 1);
title('1. Raw Received (Phase Rotation)');
xlabel('In-Phase'); ylabel('Quadrature');
axis([-2.5 2.5 -2.5 2.5]); axis square; grid on;

% Plot 2: ZF
subplot(1, 3, 2);
plot(Rx_ZF_Vis(:), '.', 'Color', [0.9 0.2 0.2], 'MarkerSize', 2);
hold on;
plot(ideal_pts, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5);
title(['2. Zero Forcing @ ' num2str(SNR_Vis) 'dB']);
subtitle('Noise Amplification');
xlabel('In-Phase'); ylabel('Quadrature');
axis([-2.5 2.5 -2.5 2.5]); axis square; grid on;

% Plot 3: MMSE
subplot(1, 3, 3);
plot(Rx_MMSE_Vis(:), '.', 'Color', [0.2 0.6 0.2], 'MarkerSize', 2);
hold on;
plot(ideal_pts, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5);
title(['3. MMSE @ ' num2str(SNR_Vis) 'dB']);
subtitle('Noise Suppression');
xlabel('In-Phase'); ylabel('Quadrature');
axis([-2.5 2.5 -2.5 2.5]); axis square; grid on;

fprintf('Simulation Complete. Figures Generated.\n');