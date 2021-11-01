% create a spectrogram/sonograph of the input signal (rect window; no-overlap)
function s = sonograph(sig, seglen, fs, fig)

N  = length(sig);

% reshape signal as column vectors of segments
sig = [sig, zeros(1, seglen-rem(length(sig), seglen))]; % zero padded signal
sig = reshape(sig, seglen, length(sig)/seglen);         % non-overlapping segment vectors

% create sonograph
s = abs(fftshift(fft(sig,[],1), 1)); % fft of each column/segment
s = s(ceil(seglen/2):seglen, :, :);  % discard negative frequencies
s = flipud(s); % flip up-down to naturally orient the frequency axis

% plot
T  = 5;  % time tick scale (number of cols per tick)
F  = 12; % frquency tick scale (number of rows per tick)
df = fs/(2*size(s,1)); % df per row
dt = seglen/fs;        % dt per column
t  = 0:T*dt:N/fs;      % time range
fr = fs/2:-F*df:0;     % frequency range

figure(fig);
imagesc(s); c = colorbar;
xticks(0.5:T:size(s,2)+0.5);
xticklabels(round(t, 2));
if length(t) > 10
    xtickangle(45);
end
yticks(0.5:F:size(s,1)+0.5);
yticklabels(round(fr/1000, 2));
if length(fr) > 10
    ytickangle(45);
end

xlabel('Time [s]', 'Interpreter', 'latex');
ylabel('Frequency [kHz]', 'Interpreter', 'latex');
ylabel(c, 'DFT Magnitude', 'Interpreter', 'latex'); 
c.Label.Position(1) = 2.5;
