% create a spectrogram/sonograph of the input signal (rect window; no-overlap)
function s = sonograph(sig, seglen)

% reshape signal as column vectors of segments
sig = [sig, zeros(1, seglen-rem(length(sig), seglen))]; % zero padded signal
sig = reshape(sig, seglen, length(sig)/seglen);         % non-overlapping segment vectors

% create sonograph
s   = abs(fftshift(fft(sig,[],1), 1)); % fft of each column/segment
s   = s(ceil(seglen/2):seglen, :, :);  % discard negative frequencies
s   = flipud(s); % flip up-down to naturally orient the frequency axis
