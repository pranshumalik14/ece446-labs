% generates a noise signal for a given order of spectral-decay
function [ns, Ns] = noisegen(n, Fs, dur) % n = 0 (white), 1 (pink), 2 (brown)

N = floor(dur*Fs); % number of sampled points of the signal (and its dft)
M = floor(N/2+1);  % number of unique fft points
f = 1:M;           % frequency indexes

Ns = exp(1i*2*pi*rand(1, M));   % initializing left half with random phase
Ns = Ns./sqrt(f.^n);            % amplitude is proportional to 1/sqrt(f^n), acc to color noise power law
Ns = [Ns conj(Ns(end-1:-1:2))]; % right half is conjugate copy of left; reconstruct the whole spectrum

% extract corresponding time-domain noise signal in [-1, 1]
ns = real(ifft(Ns));
ns = rescale(ns, -1, 1);

Ns = fftshift(fft(ns)); % new fft after rescaling
