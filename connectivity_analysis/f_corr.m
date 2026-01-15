%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              f_corr
% author - Brad Rauscher (created 2023)
% 
% Calculates correlation (pearson's coefficient) between matrices sig1 and 
% sig2 along 'dim' dimension. Removes average before calculating 
% correlation. Option to plot result as a correlation map.
% 
% INPUTS: f_corr(sig1, sig2, dim, _)
%   sig1: first signal
%   sig2: second signal
%   dim: dimension to calculate correlation on
% 
% OPTIONAL INPUTS:
%   plot: plot results as a correlation map (default = false)
% 
% OUTPUTS:
%   r: correlation value(s)
%   f: figure handle
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [r, f] = f_corr(sig1, sig2, dim, varargin)

% parse inputs
p = inputParser;
addParameter(p, 'plot', false); % plot results as correlation map

parse(p, varargin{:});

% remove mean from sig1 and sig2
sig1 = sig1 - mean(sig1, dim);
sig2 = sig2 - mean(sig2, dim);

% calculate pearson's coefficient along 'dim'
std1 = sum(sig1.^2, dim);
std2 = sum(sig2.^2, dim);
std = sqrt(std1 .* std2);
cov = sum((sig1 .* sig2), dim);

r = cov ./ std; 

% plot results
if p.Results.plot
    f = figure;
    imagesc(r, AlphaData =~ isnan(r));
    axis image off;
    c = colorbar;
    c.Label.String = 'r';
    set(gca, FontSize=14);
else
    f = [];
end