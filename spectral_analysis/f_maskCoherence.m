%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           f_maskCoherence
% author - Brad Rauscher (created 2025)
% 
% Calculates coherence between each pixel in 'sig1' and 'sig2' using the
% coherencyc function in the Chronux toolbox. Averages the calculated
% coherence values within each mask in 'masks'. 
% 
% INPUTS: f_maskCoherence(sig1, sig2, fs, masks, _)
%   sig1: first signal (H x W x T)
%   sig2: second signal  (H x W x T)
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
%   C: coherence value(s)
%   phi: phase value(s) for C
%   f: frequency value(s) for C
%   f_handle: figure handle (only outputs when plot = true)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [C, phi, f, f_handle] = ...
    f_maskCoherence(sig1, sig2, fs, masks, varargin)

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

dim = size(sig1);

N = size(masks, 3);
masks = reshape(masks, dim(1) * dim(2), N);
masks = masks(mask_all(:), :);

sig1 = reshape(sig1, dim(1) * dim(2), dim(3));
sig1 = sig1(mask_all(:), :)';

sig2 = reshape(sig2, dim(1) * dim(2), dim(3));
sig2 = sig2(mask_all(:), :)';

% calculate coherence of each column in 'sig1' and 'sig2'
params = struct;
params.Fs = fs;
params.tapers = p.Results.tapers;
params.trialave = 0;

[C_all, phi_all, ~, ~, ~, f] = coherencyc(sig1, sig2, params);
f = f';

% calculate mean C and phi for each mask
C = zeros(numel(f), N);
phi = C;

for i = 1 : N
    C(:, i) = mean(C_all(:, masks(:, i)), 2);
    phi(:, i) = mean(phi_all(:, masks(:, i)), 2);
end

% plot C and phi
if p.Results.plot
    f_handle = figure(Position = [100, 100, 700, 600]);
    tiledlayout(2, 1, ...
        TileSpacing = 'compact', ...
        Padding = 'compact');

    ax(1) = nexttile;
    imagesc(C', ...
        XData = f);

    c = colorbar(XTick = 0 : 0.1 : 1, ...
        XTickLabel = {0, '', '', '', '', '', '', '', '', '', 1});
    c.Label.String = 'coherence';

    set(ax(1), ...
        YTick = 1 : N, ...
        YTickLabel = p.Results.labels, ...
        XTick = [], ...
        FontSize = 14);

    clim([0, 1]);
    xlim(p.Results.fRange);

    ax(2) = nexttile;
    imagesc(phi', ...
        XData = f);
    c = colorbar(XTick = pi * [-1, 0, 1], XTickLabel = {'-\pi', 0, '\pi'});
    c.Label.String = 'phase (rad)';
    
    xlabel('Frequency (Hz)');
    set(ax(2), ...
        YTick = 1 : N, ...
        YTickLabel = p.Results.labels, ...
        XTick = linspace(p.Results.fRange(1), p.Results.fRange(2), 11), ...
        FontSize = 14);
    colormap cmpbbr;
    clim(pi * [-1, 1]);

    xlim(p.Results.fRange);

    colormap(ax(1), cmpinf);
else
    f_handle = [];
end