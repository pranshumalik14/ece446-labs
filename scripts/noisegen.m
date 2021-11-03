% generates a noise signal for a given order of spectral-decay
function [ns, Ns] = noisegen(c, n, Fs, dur) % n = 0 (white), 1 (pink), 2 (brown)

N = floor(dur*Fs); % number of sampled points of the signal (and its dft)
M = floor(N/2+1);

Ns = c*exp(1i*2*pi*rand(1, M));
% prepare a vector with frequency indexes 
NumUniquePts = M;     % number of the unique fft points
k = 1:NumUniquePts;         % vector with frequency indexes 
% manipulate the left half of the spectrum so the PSD
% is proportional to the frequency by a factor of 1/f, 
% i.e. the amplitudes are proportional to 1/sqrt(f)
Ns = Ns(1:NumUniquePts);      
Ns = Ns./sqrt(k.^n);
% prepare the right half of the spectrum - a conjugate copy of the left
% one except the DC component and the Nyquist component - they are unique,
% and reconstruct the whole spectrum
Ns = [Ns conj(Ns(end-1:-1:2))];
% IFFT
Ns = Ns(1:N);
ns = real(ifft(Ns));

% form the noise matrix and ensure unity standard 
% deviation and zero mean value (columnwise)
ns = bsxfun(@minus, ns, mean(ns));
ns = bsxfun(@rdivide, ns, std(ns));