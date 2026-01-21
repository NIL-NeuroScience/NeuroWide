%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           f_maskHemLag
% author - Brad Rauscher (created 2025)
% 
% Calculates cross-correlation between pixels of matrices sig1 and 
% sig2 along third dimension and averages within each mask in 'masks'. 
% 
% INPUTS: f_maskHemLag(sig1, sig2, maxlag, masks, _)
%   sig1: first signal (H x W x T)
%   sig2: second signal (H x W x T)
%   maxlag: max negative and positive lag
%   masks: masks to average within (H x W x N_masks)
% 
% OPTIONAL INPUTS:
%   detrend: remove linear trend for each pixel (default = 1)
%   plot: plots results (default = false)
%   ylabel: labels for 'masks' (cell array)
% 
% OUTPUTS:
%   r: correlation value(s)
%   lag: lag value(s) for r
%   f_handle: figure handle (only outputs when plot = true)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [r, lag, f_handle] = ...
    f_maskHemLag(sig1, sig2, maxlag, masks, varargin)

% parse inputs
p = inputParser;
addParameter(p, 'detrend', true); % plot results as correlation map
addParameter(p, 'plot', false); % plot results as correlation map
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

% calculate cross-correlation of each column in 'sig1' and 'sig2'

r_all = f_xcorr(sig1, sig2, maxlag);
lag = -maxlag : maxlag;

% calculate mean r for each mask
r = zeros(numel(lag), N);

for i = 1 : N
    r(:, i) = mean(r_all(:, masks(:, i)), 2);
end

% plot r and average r
if p.Results.plot
    f_handle = figure(Position = [100, 100, 700, 600]);
    tiledlayout(2, 1, ...
        TileSpacing = 'compact', ...
        Padding = 'compact');

    ax(1) = nexttile;
    imagesc(r', ...
        XData = lag);

    c = colorbar(XTick = -1 : 0.2 : 1, ...
        XTickLabel = {-1, '', '', '', '', '', '', '', '', '', 1});
    c.Label.String = 'r';
    colormap cmpbbr;

    set(ax(1), ...
        YTick = 1 : N, ...
        YTickLabel = p.Results.labels, ...
        XTick = [], ...
        FontSize = 14);

    clim([-1, 1]);
    xlim(maxlag * [-1, 1]);

    ax(2) = nexttile;
    plot(lag, mean(r, 2));
    
    xlabel('lag (frames)');
    set(ax(2), ...
        YTick = -1 : 1, ...
        XTick = linspace(lag(1), lag(end), 11), ...
        FontSize = 14);
    
    xlim(maxlag * [-1, 1]);
    ylim([-1 ,1]);
    ylabel('r');
    box off;

else
    f_handle = [];
end