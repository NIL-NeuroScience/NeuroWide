%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             f_maskSpectra
% author - Brad Rauscher (created 2025)
% 
% Calculates power spectral density for each pixel in 'sig' using the
% mtspectrumc function in the Chronux toolbox. Averages the calculated
% power spectral desnity values within each mask in 'masks'.
% 
% INPUTS: f_maskSpectra(sig, fs, masks, _)
%   sig: signal (H x W x T)
%   fs: sampling rate of signals (Hz)
%   masks: masks to average within (H x W x N_masks)
% 
% OPTIONAL INPUTS:
%   tapers: multitaper parameters [time-bandwidth product, N tapers]
%       (default = [5, 9])
%   plot: plots results (default = false)
%   fRange: frequency range to plot
%   ylabel: labels for 'masks' (cell array)
% 
% OUTPUTS:
%   S: spectral power density value(s)
%   f: frequency value(s) for S
%   f_handle: figure handle (only outputs when plot = true)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [S, f, f_handle] = f_maskSpectra(sig, fs, masks, varargin)

% parse inputs
p = inputParser;
addParameter(p, 'tapers', [5, 9]); % multitaper parameters
addParameter(p, 'plot', false); % plot results as correlation map
addParameter(p, 'fRange', [0, fs / 2]); % frequency range to plot
addParameter(p, 'labels', {}); % labels for 'masks'

parse(p, varargin{:});

% process inlputs
masks(isnan(masks)) = 0;
masks = logical(masks);

mask_all = logical(sum(masks, 3));

dim = size(sig);

N = size(masks, 3);
masks = reshape(masks, dim(1) * dim(2), N);
masks = masks(mask_all(:), :);

sig = reshape(sig, dim(1) * dim(2), dim(3));
sig = sig(mask_all(:), :)';

% calculate coherence of each column in 'sig1' and 'sig2'
params = struct;
params.Fs = fs;
params.tapers = p.Results.tapers;
params.trialave = 0;

[S_all, f] = mtspectrumc(sig, params);
f = f';

% calculate mean C and phi for each mask
S = zeros(numel(f), N);

for i = 1 : N
    S(:, i) = mean(S_all(:, masks(:, i)), 2);
end

% plot C and phi
if p.Results.plot
    f_handle = figure(Position = [100, 100, 700, 600]);
    tiledlayout(2, 1, ...
        TileSpacing = 'compact', ...
        Padding = 'compact');

    S_plot = log10(S);
    ax(1) = nexttile;
    imagesc(S_plot', ...
        XData = f);
    c_bounds = prctile(S_plot(:),[0.1, 99.9]);
    c_bounds_r = round(c_bounds, 2);

    c = colorbar(XTick = linspace(c_bounds(1), c_bounds(2), 5), ...
        XTickLabel = {c_bounds_r(1), '', '', '', c_bounds_r(2)});
    c.Label.String = 'log_1_0(PSD)';
    colormap cmpinf;

    set(ax(1), ...
        YTick = 1 : N, ...
        YTickLabel = p.Results.labels, ...
        XTick = [], ...
        FontSize = 14);

    clim(c_bounds);
    xlim(p.Results.fRange);

    ax(2) = nexttile;
    plot(f, mean(S_plot, 2));
    
    xlabel('f (Hz)');
    set(ax(2), ...
        YTick = linspace(c_bounds(1), c_bounds(2), 5), ...
        YTickLabel = {c_bounds_r(1), '', '', '', c_bounds_r(2)}, ...
        XTick = linspace(p.Results.fRange(1), p.Results.fRange(2), 11), ...
        FontSize = 14);
    
    xlim(p.Results.fRange);
    ylim(c_bounds);
    ylabel('log_1_0(PSD)');
    box off;
else
    f_handle = [];
end