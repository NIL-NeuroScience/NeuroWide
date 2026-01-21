%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              f_morlet
% author - Brad Rauscher (created 2024)
% 
% Spectral analysis function to calculate the instantaneous power spectral 
% density of the input signal using the Morlet wavelet approach.
% 
% INPUTS: f_morlet(y, fs, fwin, fsteps, _)
%   y: signal to estimate power of
%   fs: sampling rate of y (Hz)
%   fwin: frequency window in the form [f_low, f_high] (Hz)
%   fsteps: number of steps to take in frequency space (100 is normal)
% 
% OPTIONAL INPUTS:
%   plot: plot spectrogram
% 
% OUTPUTS:
%   spg: power spectral density spectrogram
%   t: time axis
%   fSpectogram: frequency axis
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [spg, t, fSpectogram] = f_morlet(y, fs, fwin, fsteps, varargin)

% handle inputs

p = inputParser;
addParameter(p, 'plot', 0);

parse(p, varargin{:});

% create time and frequency axes
t = 1 / fs : 1 / fs : size(y, 1) / fs;

fSpectogram=logspace(log10(fwin(1)), log10(fwin(2)), fsteps);
morletFWHM=5;

%
Ts = 1 / fs;

sigma_tc = morletFWHM / sqrt(8 * log(2));
sigma_t = sigma_tc ./ fSpectogram;
nscales = length(fSpectogram);
precision = 3;

nx = size(y, 1);
ntimes = size(y, 2);
P = zeros(nx, nscales, ntimes);

for s = 1 : nscales
    xval = (-precision * sigma_t(s) : 1 / fs : precision * sigma_t(s))';
    W = sqrt(fSpectogram(s)) * morlet_wavelet(fSpectogram(s) * xval, ...
        sigma_tc);
    P(:, s, :) = conv2(y, W, 'same') * Ts; 
end

spg = abs(P).^2;

% plot spg

if ~p.Results.plot
    return
end

% find ticks
y_int = log10(fSpectogram(1));
M = (fsteps - 1) / (log10(fSpectogram(end)) - y_int);

fun = @(x) 1 + M * (log10(x) - y_int);

[major_ticks, minor_ticks] = log_ticks(fwin);

% plot spg
plot_spg = log10(flipud((spg)'));

ax = gca;
imagesc(plot_spg, YData=size(spg, 2) : -1 : 1, ...
    XData=(1 : size(spg, 1)) / fs);
clim(prctile(plot_spg(:), [0.1, 99.9]));

c = colorbar;
c.Label.String = 'Log_1_0(Power)';
set(ax, ...
    YDir = 'normal', ...
    YTick = fun(major_ticks), ...
    YTickLabel = major_ticks, ...
    TickDir = 'out');
ax.YAxis.MinorTickValues = fun(minor_ticks);
ax.YMinorTick = 'on';

xlabel('Time (s)');
ylabel('F (Hz)');

box off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EXTRA FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create morlet wavelet for a given sigma value
function W = morlet_wavelet(t, sigma_tc)
    W = (sigma_tc * sqrt(pi))^(-0.5) * exp(-(t.^2) ...
        / (2 * sigma_tc^2)) .* exp(1i * 2 * pi * t);
end

% handle log frequency ticks for plotting
function [major_ticks, minor_ticks] = log_ticks(fwin)
    % Frequency range
    fmin = fwin(1);
    fmax = fwin(2);
    
    % Determine the decades covering the range
    decades = floor(log10(fmin)) : ceil(log10(fmax));
    
    % Initialize arrays
    major_ticks = [];
    minor_ticks = [];
    
    for d = decades
        % Minor ticks: 1-9 * 10^d
        candidates = (1 : 9) * 10^d;
        
        % Only keep within range
        candidates = candidates(candidates >= fmin & candidates <= fmax);
        
        % Add to minor ticks
        minor_ticks = [minor_ticks, candidates];
        
        % Major tick: exact powers of 10
        if 10^d >= fmin && 10^d <= fmax
            major_ticks = [major_ticks, 10^d];
        end
    end
    
    % Remove major ticks from minor ticks
    minor_ticks = setdiff(minor_ticks, major_ticks);
end

end