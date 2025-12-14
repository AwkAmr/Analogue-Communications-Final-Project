%% Analog Communication - Experiment 3: NBFM (Narrowband Frequency Modulation)
% This script implements narrowband FM modulation and demodulation
clear all;
close all;
clc;

%% Step 1: Read and Process Audio Signal
% Read the audio file
[m_t, Fs_original] = audioread('eric.wav');
m_t = m_t(:,1); % Take only one channel if stereo
t_original = (0:length(m_t)-1)/Fs_original; % Time vector

% Plot original signal in time domain
figure('Name', 'Original Audio Signal');
subplot(2,1,1);
plot(t_original, m_t); % Plot amplitude vs time
title('Original Audio Signal - Time Domain');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

% Compute spectrum of original signal
M_f = fftshift(fft(m_t)); % FFT and shift zero frequency to center
freq_original = linspace(-Fs_original/2, Fs_original/2, length(M_f)); % Frequency axis

subplot(2,1,2);
plot(freq_original/1000, abs(M_f)); % Plot magnitude spectrum in kHz
title('Original Audio Signal - Frequency Domain');
xlabel('Frequency (kHz)');
ylabel('Magnitude');
grid on;
xlim([-10 10]); % Focus on main frequency components

%% Step 2: Band-limit the signal to 4 kHz using ideal filter
% Create ideal low-pass filter (frequency domain)
cutoff_freq = 4000; % 4 kHz cutoff
ideal_filter = (abs(freq_original) <= cutoff_freq)'; % Logical mask

% Apply filter in frequency domain (element-wise multiplication)
M_f_filtered = M_f .* ideal_filter;

% Plot filtered spectrum
figure('Name', 'Filtered Signal');
subplot(2,1,1);
plot(freq_original/1000, abs(M_f_filtered));
title('Filtered Signal Spectrum (BW = 4 kHz)');
xlabel('Frequency (kHz)');
ylabel('Magnitude');
grid on;
xlim([-10 10]);

%% Step 3: Convert back to time domain
m_t_filtered = real(ifft(ifftshift(M_f_filtered))); % Inverse FFT and take real part

subplot(2,1,2);
plot(t_original, m_t_filtered);
title('Filtered Signal - Time Domain');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

%% Step 4: Play the filtered audio
fprintf('Playing filtered audio...\n');
sound(m_t_filtered, Fs_original); % Play audio at original sampling rate
pause(length(m_t_filtered)/Fs_original + 1); % Wait until audio finishes

%% Step 5: Resample for modulation
% Set carrier frequency and new sampling frequency
fc = 100e3; % Carrier frequency 100 kHz
Fs_new = 5 * fc; % New sampling rate 5*fc = 500 kHz to satisfy Nyquist

% Resample the filtered signal to new sampling rate
m_t_resampled = resample(m_t_filtered, Fs_new, Fs_original);
t_new = (0:length(m_t_resampled)-1)/Fs_new; % New time vector

fprintf('Original Fs: %d Hz\n', Fs_original);
fprintf('New Fs: %d Hz\n', Fs_new);
fprintf('Carrier Fc: %d Hz\n', fc);

%% Step 6: Generate NBFM Signal
% NBFM Theory:
% For NBFM, the modulation index β must be << 1 (typically β < 0.3)
% s_FM(t) = Ac*cos(2πfc*t + 2π*kf*∫m(τ)dτ)
% For NBFM (β << 1): s_NBFM(t) ≈ Ac*cos(2πfc*t) - Ac*β*m(t)*sin(2πfc*t)

% Normalize message signal to [-1,1]
m_t_norm = m_t_resampled / max(abs(m_t_resampled));

% Set frequency deviation (choose small value for NBFM)
% For NBFM: β = Δf / fm << 1
% Let's choose Δf = 750 Hz, fm = 4 kHz (max message frequency)
% This gives β = 0.1875 which is well within NBFM range
delta_f = 750; % Frequency deviation in Hz
kf = 2 * pi * delta_f; % Frequency sensitivity constant

% Calculate modulation index
beta = delta_f / cutoff_freq; % β = Δf / fm
fprintf('Modulation Index β = %.4f\n', beta);
fprintf('For NBFM, β should be < 0.3. Current β = %.4f\n', beta);

if beta >= 0.3
    warning('β >= 0.3! This is approaching wideband FM.');
end

% Generate carrier signal
carrier = cos(2*pi*fc*t_new');

% Generate NBFM signal using INTEGRAL METHOD (proper FM generation)
% Integrate the message signal: ∫m(τ)dτ
m_integral = cumsum(m_t_norm) / Fs_new;

% FM signal: s(t) = Ac*cos(2πfc*t + 2π*kf*∫m(τ)dτ)
Ac = 1; % Carrier amplitude
s_nbfm = Ac * cos(2*pi*fc*t_new' + kf * m_integral); % FM signal

% Plot NBFM signal in time domain (first 1000 samples)
figure('Name', 'NBFM Signal - Time Domain');
plot(t_new(1:1000), s_nbfm(1:1000));
title('NBFM Modulated Signal - Time Domain (First 1000 samples)');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

%% Step 7: Plot NBFM Spectrum
S_nbfm = fftshift(fft(s_nbfm)); % FFT and shift zero frequency
freq_new = linspace(-Fs_new/2, Fs_new/2, length(S_nbfm)); % Frequency axis

figure('Name', 'NBFM Signal - Frequency Domain');
plot(freq_new/1000, abs(S_nbfm));
title('NBFM Signal Spectrum');
xlabel('Frequency (kHz)');
ylabel('Magnitude');
grid on;
xlim([95 105]); % Focus on carrier region (100 kHz ± 5 kHz)

% Answer to Question 2: What can you make out of the resulting plot?
fprintf('\n--- Observation of NBFM Spectrum ---\n');
fprintf('The NBFM spectrum shows:\n');
fprintf('1. A dominant carrier component at fc = 100 kHz\n');
fprintf('2. Sidebands around the carrier (similar to DSB-SC pattern)\n');
fprintf('3. Multiple frequency components due to FM modulation\n');
fprintf('4. For NBFM with small β, spectrum resembles AM/DSB\n');
fprintf('5. Approximate bandwidth ≈ 2*(Δf + fm) = 2*(%.1f + 4) = %.1f kHz (Carson''s Rule)\n', delta_f/1000, 2*(delta_f+cutoff_freq)/1000);
fprintf('6. The "hill" shape shows the distribution of energy across frequencies\n');

%% Step 8: Answer to Question 3 - Condition for NBFM
fprintf('\n--- Condition for NBFM ---\n');
fprintf('The condition needed to achieve NBFM is:\n');
fprintf('β = Δf/fm << 1 (typically β < 0.3)\n');
fprintf('Where:\n');
fprintf('  Δf = frequency deviation = kf*max(|m(t)|)/(2π)\n');
fprintf('  fm = maximum message frequency\n');
if beta < 0.3
    fprintf('Current β = %.4f ✓ (NBFM satisfied)\n', beta);
else
    fprintf('Current β = %.4f ✗ (Not NBFM)\n', beta);
end

%% Step 9: Demodulate NBFM using Differentiator + Envelope Detector

% Step 9a: Differentiate the NBFM signal
% For FM: s(t) = Ac*cos(2πfc*t + 2π*kf*∫m(τ)dτ)
% Differentiation gives: ds/dt ∝ instantaneous frequency ∝ m(t)
s_diff = diff(s_nbfm) * Fs_new; % Multiply by Fs to scale properly

% Pad to maintain length
s_diff = [s_diff; s_diff(end)];

% Plot differentiated signal
figure('Name', 'Demodulation Process');
subplot(3,1,1);
plot(t_new(1:5000), s_diff(1:5000));
title('After Differentiation');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

% Step 9b: Apply Envelope Detector
% The envelope detector extracts the envelope of the differentiated signal
envelope = abs(hilbert(s_diff));

subplot(3,1,2);
plot(t_new(1:5000), envelope(1:5000));
title('After Envelope Detection');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

% Step 9c: Remove DC component and low-pass filter
% Remove DC (carrier frequency component)
dc_level = mean(envelope);
m_demod = envelope - dc_level;

% Low-pass filter to remove high-frequency components
% Design ideal LPF at message bandwidth (4 kHz)
M_demod_freq = fftshift(fft(m_demod));
lpf = (abs(freq_new) <= cutoff_freq)'; % Ideal LPF
M_demod_filtered = M_demod_freq .* lpf;
m_demod_filtered = real(ifft(ifftshift(M_demod_filtered)));

subplot(3,1,3);
plot(t_new(1:5000), m_demod_filtered(1:5000));
title('Demodulated Signal (After LPF)');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

% Normalize demodulated signal
m_demod_norm = m_demod_filtered / max(abs(m_demod_filtered));

% Downsample back to original sampling frequency for playback
m_received = resample(m_demod_norm, Fs_original, Fs_new);

% Ensure same length as original
if length(m_received) > length(m_t_filtered)
    m_received = m_received(1:length(m_t_filtered));
else
    m_received = [m_received; zeros(length(m_t_filtered)-length(m_received), 1)];
end

%% Step 10: Compare original and demodulated signals
figure('Name', 'Comparison: Original vs Demodulated');
subplot(2,1,1);
plot(t_original, m_t_filtered/max(abs(m_t_filtered)));
hold on;
plot(t_original, m_received);
title('Original (blue) vs Demodulated (orange) Signal');
xlabel('Time (s)');
ylabel('Normalized Amplitude');
legend('Original', 'Demodulated');
grid on;

% Compute spectra for comparison
M_orig_freq = fftshift(fft(m_t_filtered));
M_recv_freq = fftshift(fft(m_received));
subplot(2,1,2);
plot(freq_original/1000, abs(M_orig_freq)/max(abs(M_orig_freq)));
hold on;
plot(freq_original/1000, abs(M_recv_freq)/max(abs(M_recv_freq)));
title('Spectrum: Original vs Demodulated');
xlabel('Frequency (kHz)');
ylabel('Normalized Magnitude');
xlim([-6 6]);
legend('Original', 'Demodulated');
grid on;

%% Step 11: Play demodulated audio
fprintf('\nPlaying demodulated audio...\n');
sound(m_received, Fs_original);

%% Summary and Conclusions
fprintf('\n========== EXPERIMENT 3 SUMMARY ==========\n');
fprintf('1. NBFM modulation requires β << 1 (β = %.4f)\n', beta);
fprintf('2. NBFM spectrum resembles DSB-SC when β is small\n');
fprintf('3. Demodulation: Differentiator + Envelope Detector\n');
fprintf('4. NBFM provides better noise immunity than AM\n');
fprintf('5. Bandwidth efficiency: NBFM BW ≈ 2Δf, similar to DSB\n');
fprintf('==========================================\n');
