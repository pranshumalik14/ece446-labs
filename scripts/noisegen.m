% generates a noise signal for a given order of spectral-decay
function [ns, Ns] = noisegen(c, n, Fs, dur) % n = 0 (white), 1 (pink), 2 (brown)

N = floor(dur*Fs); % number of sampled points of the signal (and its dft)
M = floor(N/2+1);

% Ns = zeros(1,N);
% % assume a sampling rate of 10000 Hz; we are creating 
% % a one second noise clip which requires 10000 samples 
% % n determine spectral characteristics of noise 
% % n=0 is white noise, n=1 is pink, n=2 is brown 
% Ns(1)=0; % this ensures zero DC offset in signal
% for k=2:M
%     theta=2*pi*rand;
%     Ns(k)=c*exp(j*theta)/k^(n/2);
% end
% for k=M+1:N
%     Ns(k)=conj(Ns(N-k+2)); % conj is conjugate
% end
% 
% ns=ifft(Ns); % sometimes better to call x=real(ifft(fx))
% % round off errors may lead to x being complex

% % generate white noise sequence
% x = randn(1, M);
% FFT
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
ns = ifft(Ns);

% form the noise matrix and ensure unity standard 
% deviation and zero mean value (columnwise)
ns = bsxfun(@minus, ns, mean(ns));
ns = bsxfun(@rdivide, ns, std(ns));