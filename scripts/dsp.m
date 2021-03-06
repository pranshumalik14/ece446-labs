%% ece446: sensory communication
%  dsp problem set
%  author: Pranshu Malik
%  date: 4 November 2021

clc;
close all;
clear all;

%%  problem 1: manipulation of signals in the frequency domain

Fs  = 44100; % sampling frequency, in hertz
dur = 3;     % duration of the signal, in seconds
N   = floor(dur*Fs);    % number of sampled points of the signal (for dft)
df  = Fs/N;             % frequency increment in nyquist range
fr  = -Fs/2:df:Fs/2-df; % frequency range (nyquist range)

f1 = 440; % frequency of first signal, in hertz
f2 = 660; % frequency of second signal, in hertz

mag        = 0.5;  % magnitude of the two sinusoids (or their phasors)
phaseshift = pi/2; % two signals are 90deg out of phase (in time)

% creating the frequency domain signal/representation (dtft)
% delta(x): approximate delta function in the discrete sense; outputs one if abs(arugment) < eps
% F(e^jw):  dtft of cos(w0*t), which is a function of normalized angular frequency w = 2pi*f/fs
% X1(e^jw): dtft x1[n] = cos(2*pi*f1*n)
% X2(e^jw): dtft x2[n] = cos(2*pi*f2*n + phaseshift); uses the time-shift property of c/dtft
% X(e^jw):  dtft x[n]  = x1[n] + x2[n]
eps   = 1e-6;
delta = @(x) abs(x) < eps;
Fejw  = @(w, w0) mag * (pi*Fs/2) * (delta((w - w0)*Fs) + delta((w + w0)*Fs)); % todo: see why we need a 1/2 in the scale!!!
Xejw1 = @(w) Fejw(w, 2*pi*f1/Fs);
Xejw2 = @(w) exp(1i*w*Fs*phaseshift/(2*pi*f2)) * Fejw(w, 2*pi*f2/Fs);
Xejw  = @(w) Xejw1(w) + Xejw2(w);

% sampling the dtft to get dft/fft
Xk = zeros(1, N);
for k = 1:N
    Xk(k) = Xejw(2*pi*(k-1)/N - pi);
end

% inverse dft to get the sequence (time domain, sampled signal)
xn = ifft(ifftshift(Xk)); assert(isreal(xn)); % xn known to be real beforehand; assert result is real

fig_1 = figure('Name', 'Frequency to Time Domain Signal Construction', 'NumberTitle', 'off');
figure(fig_1);
lim = floor(N/64);
stem(fr(N/2-lim:N/2+lim), abs(Xk(N/2-lim:N/2+lim)));
xlabel('Frequency [Hz]', 'Interpreter', 'latex');
ylabel('DFT Magnitude, $|X[k]|$', 'Interpreter', 'latex');

fig_2 = figure('Name', 'Time Domain Signal', 'NumberTitle', 'off');
figure(fig_2);
stem(1:200, xn(1:200));
xlabel('$n$', 'Interpreter', 'latex');
ylabel('Signal, $x[n]$', 'Interpreter', 'latex');
ylim([-1.1 1.1])

%% problem 2: effects of windowing signal in time domain on the frequency spectra

Fs  = 44100; % sampling frequency, in hertz
f   = 1000;  % frequency of signal, in hertz
dur = 10;    % duration of the signal, in seconds
N   = floor(dur*Fs);    % number of sampled points of the signal (for dft)
df  = Fs/N;             % frequency increment in nyquist range
fr  = -Fs/2:df:Fs/2-df; % frequency range (nyquist range)

t   = linspace(0, dur, N); % time range
x1  = sin(2*pi*f*t);       % original signal (high sampling duration)

st = 1; ter = 9;               % start and end of the zeroing-window 
x2 = x1; x2(st*Fs:ter*Fs) = 0; % cut signal (low sampling duration)

% get ffts
X1 = fftshift(fft(x1));
X2 = fftshift(fft(x2));

fig_3 = figure('Name', 'Full Duration Signal Spectra', 'NumberTitle', 'off');
fig_4 = figure('Name', 'Windowed Signal Spectra', 'NumberTitle', 'off');
fl = 960;
fh = 1040;

figure(fig_3);
plot(fr(N/2+fl*N/Fs:N/2+fh*N/Fs), abs(X1(N/2+fl*N/Fs:N/2+fh*N/Fs)));
xlim([fl fh]);
xlabel('Frequency [Hz]', 'Interpreter', 'latex');
ylabel('DFT Magnitude, $|X_1[k]|$', 'Interpreter', 'latex');

figure(fig_4);
plot(fr(N/2+fl*N/Fs:N/2+fh*N/Fs), abs(X2(N/2+fl*N/Fs:N/2+fh*N/Fs)));
xlim([fl fh]);
xlabel('Frequency [Hz]', 'Interpreter', 'latex');
ylabel('DFT Magnitude, $|X_2[k]|$', 'Interpreter', 'latex');

%% problem 3: telephone dial tones

Fs  = 10000; % sampling frequency, in hertz
dur = 1;     % duration of the signal, in seconds
N   = floor(dur*Fs);  % number of sampled points of the signal
t   = linspace(0, dur, N); % time range

fl1 = 697;   % first lower band frequency, in hertz
fl2 = 770;   % second lower band frequency, in hertz
fu1 = 1209;  % first upper band frequency, in hertz
fu2 = 1336;  % second upper band frequency, in hertz

dialtone_1 = 0.5*(cos(2*pi*fl1*t) + cos(2*pi*fu1*t));
dialtone_2 = 0.5*(cos(2*pi*fl1*t) + cos(2*pi*fu2*t));
dialtone_4 = 0.5*(cos(2*pi*fl2*t) + cos(2*pi*fu1*t));

%% problem 4: dialtone decode algorithm

% simple dtmf decode algorithm for a subset of tones (can be extended)
input_dtmf_tone = '';
[signal, Fs]    = audioread(input_dtmf_tone);

N   = length(signal);   % number of points in input signal and its fft
df  = Fs/N;             % frequency increment in nyquist range
fr  = -Fs/2:df:Fs/2-df; % frequency range (nyquist range)

fr_eps      = 10;         % frequency window around pure tones
dtmf_matrix = [1 2; 4 5]; % (partial) dtmf decode matrix; indexed by lower and upper band freqs
signal_fft  = fftshift(fft(signal)); % fft of the input signal

% get magnitude of each tone's contribution to the signal spectrum
m_fl1 = sum(2*abs(signal_fft(abs(fr - fl1) < fr_eps)));
m_fl2 = sum(2*abs(signal_fft(abs(fr - fl2) < fr_eps)));
m_fu1 = sum(2*abs(signal_fft(abs(fr - fu1) < fr_eps)));
m_fu2 = sum(2*abs(signal_fft(abs(fr - fu2) < fr_eps)));

% get dtmf decoding matrix indexes
[~, fl] = max([m_fl1 m_fl2]);
[~, fh] = max([m_fu1 m_fu2]);

dtmf_tone = dtmf_matrix(fl, fh); % algorithm output

%% problem 5: telephone event tones

Fs    = 10000; % sampling frequency, in hertz
dbus  = 10;    % duration of busy signal, in seconds
dcont = 5;     % duration of continuous dialtone, in seconds
Nbus  = floor(dbus*Fs);  % number of sampled points in busy signal
Ncont = floor(dcont*Fs); % number of sampled points in continuous dialtone
tbus  = linspace(0, dbus, Nbus);   % time range for busy signal
tcont = linspace(0, dcont, Ncont); % time range for continuous dialtone

busy_signal   = 0.5*(cos(2*pi*480*tbus) + cos(2*pi*620*tbus));
busy_signal   = busy_signal .* pulstran(tbus, 0.25:1:10, @rectpuls, 0.5);
cont_dialtone = 0.5*(cos(2*pi*350*tcont) + cos(2*pi*440*tcont));

%% problem 6: helmholtz resonator

[empty_bottle, Fs_empty] = audioread('../audio/6/empty_bottle_resonance.m4a');
[fill_bottle, Fs_fill]   = audioread('../audio/6/one_third_filled_bottle_resonance.m4a');

N_empty = length(empty_bottle); % number of points in empty-bottle resonance
N_fill  = length(fill_bottle);  % number of points in (partially) filled-bottle resonance

fft_empty  = fftshift(fft(empty_bottle)); % fft of empty-bottle resonance recording
fft_fill   = fftshift(fft(fill_bottle));  % fft of filled-bottle resonance recording

df_empty = Fs_empty/N_empty; % frequency increment in empty-bottle recording
df_fill  = Fs_fill/N_fill;   % frequency increment in filled-bottle recording
fr_empty = -Fs_empty/2:df_empty:Fs_empty/2-df_empty; % empty-bottle recording frequency range
fr_fill  = -Fs_fill/2:df_fill:Fs_fill/2-df_fill;     % filled-bottle recording frequency range
fsptr_empty = 100 < fr_empty & fr_empty < 300;       % plotting spectrum for empty bottle
fsptr_fill  = 100 < fr_fill  & fr_fill  < 300;       % plotting spectrum for filled bottle

fig_empty = figure('Name', 'Empty Bottle Resonance Spectra', 'NumberTitle', 'off');
fig_fill  = figure('Name', 'Partially Filled Bottle Resonance Spectra', 'NumberTitle', 'off');

figure(fig_empty);
plot(fr_empty(fsptr_empty), abs(fft_empty(fsptr_empty)));
xlabel('Frequency [Hz]', 'Interpreter', 'latex');
ylabel('DFT Magnitude', 'Interpreter', 'latex');

figure(fig_fill);
plot(fr_fill(fsptr_fill), abs(fft_fill(fsptr_fill)));
xlabel('Frequency [Hz]', 'Interpreter', 'latex');
ylabel('DFT Magnitude', 'Interpreter', 'latex');

%% problem 7: voice and instrument timbre comparison

Fs  = 44100;         % sampling frequency, in hertz
dur = 5;             % duration of reference tone, in seconds
N   = floor(dur*Fs); % number of points in reference tone and its fft
t   = linspace(0, dur, N);      % time range for reference tone
df  = Fs/N;                     % frequency increment in nyquist range
fr  = -Fs/2:df:Fs/2-df;         % frequency range (nyquist range)
a440_reftone = cos(2*pi*440*t); % 440Hz (middle A) reference tone

% violin and voice recordings at 440Hz
[a440_violin, Fs_violin] = audioread('../audio/7/a440_violin.m4a');
[a440_voice, Fs_voice]   = audioread('../audio/7/a440_voice.wav');

% fft of all signals
N_violin = length(a440_violin); % number of points in violin signal
N_voice  = length(a440_voice);  % number of points in voice signal
df_violin = Fs_violin/N_violin; % frequency increment in violin recording
df_voice  = Fs_voice/N_voice;   % frequency increment in voice recording
fr_violin = -Fs_violin/2:df_violin:Fs_violin/2-df_violin; % violin frequency range
fr_voice  = -Fs_voice/2:df_voice:Fs_voice/2-df_voice;     % voice frequency range

fft_reftone = fftshift(fft(a440_reftone));
fft_violin  = fftshift(fft(a440_violin));
fft_voice   = fftshift(fft(a440_voice));

% plot
fig_5 = figure('Name', 'A-440Hz Pure Tone Spectra', 'NumberTitle', 'off');
fig_6 = figure('Name', 'A-440Hz Violin Spectra', 'NumberTitle', 'off');
fig_7 = figure('Name', 'A-440Hz Voice Spectra', 'NumberTitle', 'off');
fig_8 = figure('Name', 'A-440Hz Violin Harmonics', 'NumberTitle', 'off');
fig_9 = figure('Name', 'A-440Hz Voice Harmonics', 'NumberTitle', 'off');
fsptr = 3000;
fhar  = 1500;

figure(fig_5);
plot(fr(abs(fr) < fsptr), abs(fft_reftone(abs(fr) < fsptr)));
xlabel('Frequency [Hz]', 'Interpreter', 'latex');
ylabel('DFT Magnitude', 'Interpreter', 'latex');
figure(fig_6);
plot(fr_violin(abs(fr_violin) < fsptr), abs(fft_violin(abs(fr_violin) < fsptr)));
xlabel('Frequency [Hz]', 'Interpreter', 'latex');
ylabel('DFT Magnitude', 'Interpreter', 'latex');
figure(fig_7);
plot(fr_voice(abs(fr_voice) < fsptr), abs(fft_voice(abs(fr_voice) < fsptr)));
xlabel('Frequency [Hz]', 'Interpreter', 'latex');
ylabel('DFT Magnitude', 'Interpreter', 'latex');

figure(fig_8);
plot(fr_violin(abs(fr_violin - fhar) < fhar), abs(fft_violin(abs(fr_violin - fhar) < fhar)));
xlabel('Frequency [Hz]', 'Interpreter', 'latex');
ylabel('DFT Magnitude', 'Interpreter', 'latex');
xticks(0:440:2*fhar);
xtickangle(45);
figure(fig_9);
plot(fr_voice(abs(fr_voice - fhar) < fhar), abs(fft_voice(abs(fr_voice - fhar) < fhar)));
xlabel('Frequency [Hz]', 'Interpreter', 'latex');
ylabel('DFT Magnitude', 'Interpreter', 'latex');
xticks(0:216:2*fhar);
xtickangle(45);

%% problem 8: linear chirp

Fs  = 16000; % sampling frequency, in hertz
dur = 1;     % duration of the signal, in seconds
N   = floor(dur*Fs);       % number of sampled points of the signal
t   = linspace(0, dur, N); % time range

f0  = 500;  % initial frequency of signal, in hertz
m   = 7000; % slope of instantaneous frequency, in hertz
phi = @(t) 1/2*m*t.^2 + f0*t; % time-dependent linear chirp signal phase

Fs_sub = 8000;                    % sub-sampling frequency, in hertz
N_sub  = floor(dur*Fs_sub);       % number of points in sub-sampled signal
t_sub  = linspace(0, dur, N_sub); % time range

% create chirps
x     = cos(2*pi*phi(t));     % linear chirp signal
x_sub = cos(2*pi*phi(t_sub)); % under/sub-sampled linear chirp signal

% ffts and plots
df     = Fs/N;         % frequency increment
df_sub = Fs_sub/N_sub; % frequency increment in subsampled signal
fr     = -Fs/2:df:Fs/2-df; % frequency range
fr_sub = -Fs_sub/2:df_sub:Fs_sub/2-df_sub; % frequency range in subsampled signal

fft_x     = fftshift(fft(x));     % fft of chirp
fft_x_sub = fftshift(fft(x_sub)); % fft of subsampled chirp

fig_10 = figure('Name', 'Linear Chirp', 'NumberTitle', 'off');
fig_11 = figure('Name', 'Subsampled Linear Chirp', 'NumberTitle', 'off');

figure(fig_10);
plot(fr, abs(fft_x));
xlabel('Frequency [Hz]', 'Interpreter', 'latex');
ylabel('DFT Magnitude', 'Interpreter', 'latex');

figure(fig_11);
plot(fr_sub, abs(fft_x_sub));
xlabel('Frequency [Hz]', 'Interpreter', 'latex');
ylabel('DFT Magnitude', 'Interpreter', 'latex');

%% problem 9: manually constructing a spectrogram (or sonograph)

fig_12 = figure('Name', 'Linear Chirp Sonograph', 'NumberTitle', 'off');
fig_13 = figure('Name', 'Subsampled Linear Chirp Sonograph', 'NumberTitle', 'off');

Xs     = sonograph(x, 256, Fs, fig_12);
Xs_sub = sonograph(x_sub, 256, Fs_sub, fig_13);

%% problem 10: doppler effect, spectrograms, and calculation of velocity

% compute and plot the spectrogram
fig_14 = figure('Name', 'Jet Takeoff Spectrogram and Doppler Shift Tracking', 'NumberTitle', 'off');
figure(fig_14);
[doppler, Fs_dop] = audioread('../audio/10/doppler.m4a');
spectrogram(doppler, 256*20, 256*15, 256*20, Fs_dop, 'yaxis');

[s,f,t,p] = spectrogram(doppler, 256*20, 256*15, 256*20, Fs_dop, 'yaxis');

% track the most significant doppler shift ("doppler chirp")
fr1 = f > 8500 & f < 9000; % first frequency range
tr1 = t > 1.4  & t < 2.3;  % first time range
fr2 = f > 7700 & f < 8800; % second frequency range
tr2 = t > 2.3  & t < 3.1;  % second time range
fr3 = f > 5550 & f < 7700; % third frequency range
tr3 = t > 3.1  & t < 4;    % third time range
fr4 = f > 5400 & f < 5650; % forth frequency range
tr4 = t > 4  & t < 4.8;    % fourth time range

m1 = medfreq(p(fr1, tr1), f(fr1)); % median frequency being tracked over fr1 and tr1
m2 = medfreq(p(fr2, tr2), f(fr2)); % median frequency being tracked over fr2 and tr2
m3 = medfreq(p(fr3, tr3), f(fr3)); % median frequency being tracked over fr3 and tr3
m4 = medfreq(p(fr4, tr4), f(fr4)); % median frequency being tracked over fr4 and tr4

% collate tracked median (f,t) pairs and smooth data (over a moving window)
m  = smooth([m1 m2 m3 m4]);
tm = [t(tr1) t(tr2) t(tr3) t(tr4)];

% plot the tracked shift/chirp
hold on;
pm = plot(tm, m/1000, 'linewidth', 4);
hold off;
pm.Color = [1.0, 0, 0, 0.1];

% apparent frequency profile fitting using non-linear regression on tracking data
N_dop = length(doppler);         % number of samples in the jet-takeoff doppler-effect sound signal
dur   = floor(N_dop/Fs_dop);     % duration of the signal, in seconds
t_dop = linspace(0, dur, N_dop); % corresponding timepoints in the doppler-effect signal
c     = 343; % assumed speed of sound, in m/s, while recording

% create model function explaining the data, or its expected relationship:
% f_app(ts) is the model function, where ts is time in seconds. f_app has parameters
% param = [f0, v, d, t0], where f0 is the stationary source freq, v is the source
% velocity, d is the closest distance by which source passes the mic, and t0 is the time offset.
% we set the initial guess at f0 = 7000Hz, v = 75m/s, d = 20m after some visual inspection.
f_app  = @(param, ts) param(1) .* c./(c + param(2).*sin(atan(param(2).*(ts - param(4))./param(3))));
param0 = [7000, 72, 20, 3.2];

% non-linear regression fit using the model
opts   = statset('TolFun', 1e-8);                       % set fit tolerance
mdl    = fitnlm(tm, m, f_app, param0, 'Options', opts); % non-linear regression fit
params = mdl.Coefficients{:, 'Estimate'};               % extract the estimated parameters

% output best-fit estimate parameters
f0 = params(1)
v  = params(2)
d  = params(3)

% plot the best-fit profile and tracking data
fig_15 = figure('Name', 'Apparent Non-linear Frequency Profile Fitting', 'NumberTitle', 'off');

figure(fig_15);
hold on;
pm = plot(tm, m, 'linewidth', 4);
pm.Color = [1.0, 0, 0, 0.25];
plot(t_dop, f_app(params, t_dop),'-b');
hold off;
xlabel('Time [s]', 'Interpreter', 'latex');
ylabel('Frequency [Hz]', 'Interpreter', 'latex');
set(legend('Doppler-shift Data', 'Best-fit $f_{app}(t)$'), 'Interpreter', 'latex');

%% problem 11: noise generation

Fs  = 10000; % sampling frequency, in hertz
dur = 5;     % duration of the noise signal, in seconds
N   = floor(dur*Fs);       % number of sampled points of the signal
t   = linspace(0, dur, N); % time range
df  = Fs/N;                % frequency increment in nyquist range
fr  = -Fs/2:df:Fs/2-df;    % frequency range (nyquist range)

[wn, Wn] = noisegen(0, Fs, dur); % white noise
[pn, Pn] = noisegen(1, Fs, dur); % pink noise
[bn, Bn] = noisegen(2, Fs, dur); % brown noise

% compute the noise power (amplitude squared) from fft
pWn = Wn.*conj(Wn);
pPn = Pn.*conj(Pn);
pBn = Bn.*conj(Bn);

% plot noise signal in time and frequency domains
fig_16 = figure('Name', 'White Noise over Time', 'NumberTitle', 'off');
fig_17 = figure('Name', 'Pink Noise over Time', 'NumberTitle', 'off');
fig_18 = figure('Name', 'Brown Noise over Time', 'NumberTitle', 'off');
fig_19 = figure('Name', 'White Noise Power Spectrum', 'NumberTitle', 'off');
fig_20 = figure('Name', 'Pink Noise Power Spectrum', 'NumberTitle', 'off');
fig_21 = figure('Name', 'Brown Noise Power Spectrum', 'NumberTitle', 'off');

figure(fig_16);
plot(t, wn);
xlabel('Time [s]', 'Interpreter', 'latex');
ylabel('Amplitude', 'Interpreter', 'latex');

figure(fig_17);
plot(t, pn);
xlabel('Time [s]', 'Interpreter', 'latex');
ylabel('Amplitude', 'Interpreter', 'latex');

figure(fig_18);
plot(t, bn);
xlabel('Time [s]', 'Interpreter', 'latex');
ylabel('Amplitude', 'Interpreter', 'latex');

figure(fig_19);
semilogx(fr(fr > 0), 10*log10(pWn(fr > 0)));
grid on;
xlabel('Frequency [Hz]', 'Interpreter', 'latex');
ylabel('Power [dB]', 'Interpreter', 'latex');

figure(fig_20);
semilogx(fr(fr > 0), 10*log10(pPn(fr > 0)));
grid on;
xlabel('Frequency [Hz]', 'Interpreter', 'latex');
ylabel('Power [dB]', 'Interpreter', 'latex');

figure(fig_21);
semilogx(fr(fr > 0), 10*log10(pBn(fr > 0)));
grid on;
xlabel('Frequency [Hz]', 'Interpreter', 'latex');
ylabel('Power [dB]', 'Interpreter', 'latex');

%% problem 12: uniformly random (white) noise signal

N = 1000;       % number of points in the noise signal
u = rand(1, N); % uniform random variable with support over (0,1)
w = 2*u - 1;    % (uniform) white noise

fig_22 = figure('Name', 'Uniform White Noise', 'NumberTitle', 'off');
figure(fig_22);
plot(1:N, w);
xlabel('Time Sequence [n]', 'Interpreter', 'latex');
ylabel('Amplitude', 'Interpreter', 'latex');

%% autoexport figures to (pdf) files
%  note: uncomment to save again

% savefig(fig_1, '../figs/problem1_fft');
% savefig(fig_2, '../figs/problem1_sig');
% savefig(fig_3, '../figs/problem2_fft_high_dur');
% savefig(fig_4, '../figs/problem2_fft_low_dur');
% savefig(fig_5, '../figs/problem7_a440_pure_spectra');
% savefig(fig_6, '../figs/problem7_a440_violin_spectra');
% savefig(fig_7, '../figs/problem7_a440_voice_spectra');
% savefig(fig_8, '../figs/problem7_a440_violin_harmonics');
% savefig(fig_9, '../figs/problem7_a440_voice_harmonics');
% savefig(fig_10, '../figs/problem8_fft_linear_chirp');
% savefig(fig_11, '../figs/problem8_fft_subsampled_linear_chirp');
% savefig(fig_12, '../figs/problem9_sonogram_linear_chirp');
% savefig(fig_13, '../figs/problem9_sonogram_subsampled_linear_chirp');
% savefig(fig_14, '../figs/problem10_jet_spectrogram_doppler_shift_tracking');
% savefig(fig_15, '../figs/problem10_f_app_best_fit_estimate');
% savefig(fig_16, '../figs/problem11_white_noise_time');
% savefig(fig_17, '../figs/problem11_pink_noise_time');
% savefig(fig_18, '../figs/problem11_brown_noise_time');
% savefig(fig_19, '../figs/problem11_white_noise_power_spectrum_db');
% savefig(fig_20, '../figs/problem11_pink_noise_power_spectrum_db');
% savefig(fig_21, '../figs/problem11_brown_noise_power_spectrum_db');
% savefig(fig_22, '../figs/problem12_uniform_white_noise_sequence');
% savefig(fig_empty, '../figs/problem6_empty_bottle_resonance_spectra');
% savefig(fig_fill, '../figs/problem6_filled_bottle_resonance_spectra');
