%% ece446: sensory communication
%  dsp problem set
%  author: Pranshu Malik
%  date: 28 October 2021

clc;
close all;
clear all;

% todo: save all audio/sound files in respective folders

%%  problem 1: manipulation of signals in the frequency domain

Fs  = 44100; % sampling frequency, in hertz
dur = 3;     % duration of the signal, in seconds
N   = floor(dur * Fs);  % number of sampled points of the signal (for dft)
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
xn = ifft(ifftshift(Xk)); % todo: this signal is not real before ifft! not sure why it is needed!!

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
N   = floor(dur * Fs);  % number of sampled points of the signal (for dft)
df  = Fs/N;             % frequency increment in nyquist range
fr  = -Fs/2:df:Fs/2-df; % frequency range (nyquist range)

t   = linspace(0, dur, dur*Fs + 1); % time range
x1  = sin(2*pi*f*t);                % original signal (high sampling duration)

s = 1; t = 9;               % start and end of the zeroing-window 
x2 = x1; x2(s*Fs:t*Fs) = 0; % cut signal (low sampling duration)

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
t   = linspace(0, dur, dur*Fs + 1); % time range

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
tbus  = linspace(0, dbus, dbus*Fs + 1);   % time range for busy signal
tcont = linspace(0, dcont, dcont*Fs + 1); % time range for continuous dialtone

busy_signal   = 0.5*(cos(2*pi*480*tbus) + cos(2*pi*620*tbus));
busy_signal   = busy_signal .* pulstran(tbus, 0.25:1:10, @rectpuls, 0.5);
cont_dialtone = 0.5*(cos(2*pi*350*tcont) + cos(2*pi*440*tcont));

%% problem 6: helmholtz resonator


%% problem 7: voice and instrument timbre comparison

Fs  = 44100; % sampling frequency, in hertz
dur = 5;     % duration of reference tone, in seconds
t   = linspace(0, dur, dur*Fs + 1); % time range for reference tone
N   = length(t);        % number of points in reference tone and its fft
df  = Fs/N;             % frequency increment in nyquist range
fr  = -Fs/2:df:Fs/2-df; % frequency range (nyquist range)
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

%% problem 8: chirp


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
