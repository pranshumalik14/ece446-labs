% create a spectrogram/sonograph of the input signal (rect window; no-overlap)
function s = sonograph(sig, seglen, fs)

N  = length(sig);

% reshape signal as column vectors of segments
sig = [sig, zeros(1, seglen-rem(length(sig), seglen))]; % zero padded signal
sig = reshape(sig, seglen, length(sig)/seglen);         % non-overlapping segment vectors

% create sonograph
s = abs(fftshift(fft(sig,[],1), 1)); % fft of each column/segment
s = s(ceil(seglen/2):seglen, :, :);  % discard negative frequencies
s = flipud(s); % flip up-down to naturally orient the frequency axis

% plot
df = fs/(2*size(s,1)); % df per row
dt = seglen/fs;        % dt per column
t  = 0:(5*dt):(N/fs);  % time range
fr = (fs/2):-(10*df):floor(df);  % frequency range

imagesc(s); colorbar;
xticks(0.5:5:size(s,2)+0.5);
xticklabels(t);
yticks(0.5:10:size(s,1)+0.5);
yticklabels(fr);
