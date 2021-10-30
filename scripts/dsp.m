%% ece446: sensory communication
%  dsp problem set
%  author: Pranshu Malik
%  date: 28 October 2021

clc;
close all;
clear all;

% todo: address all todos!

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

%% autoexport figures to (pdf) files
%  note: uncomment to save again

% savefig(fig_1, '../figs/problem1_fft');
% savefig(fig_2, '../figs/problem1_sig');
% savefig(fig_3, '../figs/problem2_fft_high_dur');
% savefig(fig_4, '../figs/problem2_fft_low_dur');
