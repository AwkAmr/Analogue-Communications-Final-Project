clc; clear; close all;
%% =======================
% 1. Read Audio Signal
% =======================
[fileName, pathName] = uigetfile('*.wav');
[signal, fs] = audioread(fullfile(pathName, fileName));
signal = mean(signal,2);              % Stereo to mono
signal = signal / max(abs(signal));   % Normalize
signal = signal(:);                   % Force column

N = length(signal);
t = (0:N-1)'/fs;
f = (-N/2:N/2-1)'*(fs/N);

SIGNAL = fftshift(fft(signal));

figure;
plot(f, abs(SIGNAL));
xlabel('Frequency (Hz)'); ylabel('Magnitude');
title('Message Spectrum');

%% =======================
% 2. Ideal Low-Pass Filter (4 kHz)
% =======================
cutoff = 4000;
H = abs(f) <= cutoff;

FILTERED_SIGNAL = SIGNAL .* H;
filtered_signal = ifft(ifftshift(FILTERED_SIGNAL),'symmetric');
filtered_signal = filtered_signal(:);

figure;
subplot(2,1,1);
plot(t, filtered_signal);
xlabel('Time (s)'); ylabel('Amplitude');
title('Filtered Signal (Time Domain)');

subplot(2,1,2);
plot(f, abs(FILTERED_SIGNAL));
xlabel('Frequency (Hz)'); ylabel('Magnitude');
title('Filtered Signal (Frequency Domain)');

sound(filtered_signal, fs);

%% =======================
% 3. Resampling
% =======================
Fc = 100e3;
Fs = 5*Fc;

signal_rs = resample(filtered_signal, Fs, fs);
signal_rs = signal_rs / max(abs(signal_rs));
signal_rs = signal_rs(:);

Nrs = length(signal_rs);
t_rs = (0:Nrs-1)'/Fs;
f_rs = (-Nrs/2:Nrs/2-1)'*(Fs/Nrs);

carrier = cos(2*pi*Fc*t_rs);

%% =======================
% 4. DSB-SC Modulation
% =======================
dsb_sc = signal_rs .* carrier;
DSB_SC_FFT = fftshift(fft(dsb_sc));

figure;
subplot(2,1,1);
plot(t_rs, dsb_sc);
xlabel('Time (s)'); ylabel('Amplitude');
title('DSB-SC Signal (Time Domain)');

subplot(2,1,2);
plot(f_rs, abs(DSB_SC_FFT));
xlabel('Frequency (Hz)'); ylabel('Magnitude');
title('DSB-SC Spectrum');

%% =======================
% 5. DSB-TC Modulation
% =======================
DC_bias = 2*max(abs(signal_rs));   % modulation index = 0.5
dsb_tc = (signal_rs + DC_bias) .* carrier;
DSB_TC_FFT = fftshift(fft(dsb_tc));

figure;
subplot(2,1,1);
plot(t_rs, dsb_tc);
xlabel('Time (s)'); ylabel('Amplitude');
title('DSB-TC Signal (Time Domain)');

subplot(2,1,2);
plot(f_rs, abs(DSB_TC_FFT));
xlabel('Frequency (Hz)'); ylabel('Magnitude');
title('DSB-TC Spectrum');

%% =======================
% 6. Envelope Detection
% =======================
env_sc = abs(hilbert(dsb_sc));
env_tc = abs(hilbert(dsb_tc));

figure;
subplot(2,1,1);
plot(t_rs, env_sc);
xlabel('Time (s)'); ylabel('Amplitude');
title('Envelope Detection (DSB-SC)');

subplot(2,1,2);
plot(t_rs, env_tc);
xlabel('Time (s)'); ylabel('Amplitude');
title('Envelope Detection (DSB-TC)');

%% =======================
% 7. Envelope Output (DSB-TC)
% =======================
rec_tc = env_tc - mean(env_tc);
rec_tc = resample(rec_tc, fs, Fs);
rec_tc = rec_tc / max(abs(rec_tc));

sound(rec_tc, fs);

% NOTE:
% Envelope detector works only for DSB-TC
% It fails for DSB-SC due to carrier suppression

%% =======================
% 8. Coherent Detection (DSB-SC)
% =======================
SNRs = [0 10 30];
[b,a] = butter(6, 5000/(Fs/2));

for k = 1:length(SNRs)
    noisy = awgn(dsb_sc, SNRs(k), 'measured');

    lo = 2*cos(2*pi*Fc*t_rs);
    mixed = noisy .* lo;

    baseband = filtfilt(b,a,mixed);
    rec = resample(baseband, fs, Fs);
    rec = rec / max(abs(rec));

    Nrec = length(rec);
    f_rec = (-Nrec/2:Nrec/2-1)'*(fs/Nrec);

    figure;
    subplot(2,1,1);
    plot((0:Nrec-1)'/fs, rec);
    xlabel('Time (s)'); ylabel('Amplitude');
    title(['Recovered DSB-SC (SNR = ',num2str(SNRs(k)),' dB)']);

    subplot(2,1,2);
    plot(f_rec, abs(fftshift(fft(rec))));
    xlabel('Frequency (Hz)'); ylabel('Magnitude');

    if ENABLE_SOUND
        sound(rec, fs);
        pause(Nrec/fs + 1);
    end
end

%% =======================
% 9. Frequency Offset (Carrier Frequency Offset)
% =======================
Fc_err = 100.1e3;
lo_freq = 2*cos(2*pi*Fc_err*t_rs);

mixed_freq = dsb_sc .* lo_freq;
rec_freq = filtfilt(b,a,mixed_freq);
rec_freq = resample(rec_freq, fs, Fs);
rec_freq = rec_freq / max(abs(rec_freq));

%% =======================
% 10. Phase Offset
% =======================
phi = 20*pi/180;  % 20 degrees
lo_phase = 2*cos(2*pi*Fc*t_rs + phi);

mixed_phase = dsb_sc .* lo_phase;
rec_phase = filtfilt(b,a,mixed_phase);
rec_phase = resample(rec_phase, fs, Fs);
rec_phase = rec_phase / max(abs(rec_phase));

%% =======================
% 11. Offset Results
% =======================
figure;
subplot(2,2,1);
plot(rec_freq);
title('Recovered Signal (Frequency Offset)');
xlabel('Samples'); ylabel('Amplitude');

subplot(2,2,2);
plot(abs(fftshift(fft(rec_freq))));
title('Spectrum (Frequency Offset)');

subplot(2,2,3);
plot(rec_phase);
title('Recovered Signal (Phase Offset)');
xlabel('Samples'); ylabel('Amplitude');

subplot(2,2,4);
plot(abs(fftshift(fft(rec_phase))));
title('Spectrum (Phase Offset)');

sound(rec_freq, fs);
pause(length(rec_freq)/fs + 1);
sound(rec_phase, fs);

