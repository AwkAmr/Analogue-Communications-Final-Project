% Step 1: Read audio file
[m, Fs] = audioread('eric.wav');
m = m(:,1); % Ensure mono
t = (0:length(m)-1)/Fs;

% Spectrum
N = length(m);
f = (-N/2:N/2-1)' * (Fs/N);  % column vector
M_f = fftshift(fft(m));

figure;
plot(f, abs(M_f));
xlabel('Frequency (Hz)');
ylabel('|M(f)|');
title('Original Audio Spectrum');
grid on;


% Step 2: Ideal LPF (4 kHz)
BW = 4000;
H = abs(f) <= BW;

M_filt_f = M_f .* H;


% Step 3: Time-domain signal
m_filt = real(ifft(ifftshift(M_filt_f)));

figure;
plot(t(1:2000), m_filt(1:2000));
xlabel('Time (s)');
ylabel('Amplitude');
title('Filtered Message Signal (Time Domain)');

disp('Step 3: Filtered baseband message');
sound(m_filt, Fs);
pause;


% Step 4: DSB-SC Modulation
Fc = 100e3;
Fs_mod = 5*Fc;

m_resampled = resample(m_filt, Fs_mod, Fs);
t_mod = (0:length(m_resampled)-1)' / Fs_mod; % column vector

carrier = cos(2*pi*Fc*t_mod);
dsb_sc = m_resampled .* carrier;

% Spectrum
N2 = length(dsb_sc);
f2 = (-N2/2:N2/2-1)' * (Fs_mod/N2); % column vector

DSB_f = fftshift(fft(dsb_sc));

figure;
plot(f2, abs(DSB_f));
xlabel('Frequency (Hz)');
ylabel('|DSB(f)|');
title('DSB-SC Spectrum');
grid on;


% Step 5: Ideal SSB (LSB)
H_ssb = (f2 <= Fc) & (f2 >= Fc-BW);
H_ssb = H_ssb | ((f2 >= -Fc) & (f2 <= -Fc+BW));

SSB_f = DSB_f .* H_ssb;
ssb = real(ifft(ifftshift(SSB_f)));

figure;
plot(f2, abs(SSB_f));
xlabel('Frequency (Hz)');
ylabel('|SSB(f)|');
title('SSB-LSB Spectrum (Ideal Filter)');
grid on;


% Step 6: Coherent Detection
rx = ssb .* (2*cos(2*pi*Fc*t_mod));

RX_f = fftshift(fft(rx));
H_rec = abs(f2) <= BW;

m_rec_f = RX_f .* H_rec;
m_rec = real(ifft(ifftshift(m_rec_f)));

figure;
plot(t_mod(1:3000), m_rec(1:3000));
xlabel('Time (s)');
ylabel('Amplitude');
title('Recovered Message (Ideal Coherent Detection)');

m_rec_audio = resample(m_rec, Fs, Fs_mod);  % resample audio to be able to play it
disp('Step 6: Recovered message (ideal coherent detection)');
sound(m_rec_audio, Fs);
pause;


% Step 7: Butterworth Filter
[b,a] = butter(4, BW/(Fs_mod/2));
m_rec_butter = filter(b,a,rx);

figure;
plot(t_mod, m_rec_butter);
xlabel('Time (s)');
ylabel('Amplitude');
title('Recovered Message (Butterworth Filter)');

m_rec_butter_audio = resample(m_rec_butter, Fs, Fs_mod);
disp('Step 7: Recovered message (Butterworth filter)');
sound(m_rec_butter_audio, Fs);
pause;


% Step 8: Noise Effect
SNRs = [0 10 30];

for i = 1:length(SNRs)
    noisy_ssb = awgn(ssb, SNRs(i), 'measured');
    
    rx_n = noisy_ssb .* (2*cos(2*pi*Fc*t_mod));
    RXN_f = fftshift(fft(rx_n));
    
    m_n_f = RXN_f .* H_rec;
    m_n = real(ifft(ifftshift(m_n_f)));
    
    figure;
    plot(t_mod(1:3000), m_n(1:3000));
    title(['Recovered Signal, SNR = ', num2str(SNRs(i)), ' dB']);
    xlabel('Time (s)');
    ylabel('Amplitude');

    m_n_audio = resample(m_n, Fs, Fs_mod);
    disp(['Step 8: Recovered message with SNR = ', num2str(SNRs(i)), ' dB']);
    sound(m_n_audio, Fs);
    pause;
end


% Step 9: SSB-TC
A = 2*max(abs(m_resampled));
ssb_tc = ssb + A*cos(2*pi*Fc*t_mod);

% Envelope Detection
env = abs(hilbert(ssb_tc));

env_f = fftshift(fft(env));
env_rec_f = env_f .* H_rec;
env_rec = real(ifft(ifftshift(env_rec_f)));

figure;
plot(t_mod, env_rec);
xlabel('Time (s)');
ylabel('Amplitude');
title('Recovered Message (Envelope Detection)');

env_audio = resample(env_rec, Fs, Fs_mod);
disp('Step 9: Envelope detected message (SSB-TC)');
sound(env_audio, Fs);
pause;
