%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              f_hemLag
% author - Brad Rauscher (created 2024)
% 
% Calculates cross-correlation between pixels of matrices sig1 and 
% sig2 along third dimension. 
% 
% INPUTS: f_hemLag(sig1, sig2, maxlag, mask)
%   sig1: first signal
%   sig2: second signal
%   maxlag: max negative and positive lag
%   brain_mask: mask of brain exposure (2D NaN image). Leave empty if no
%       mask is needed
% 
% OPTIONAL INPUTS:
%   detrend: remove linear trend for each pixel (default = 1)
% 
% OUTPUTS:
%   r: correlation value(s) for each lag in 'lag'
%   lag: lag values for each frame in 'xcorr'
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [r, lag] = f_hemLag(sig1, sig2, maxlag, brain_mask, varargin)

% parse inputs
p = inputParser;
addParameter(p, 'detrend', true); % plot results as correlation map

parse(p, varargin{:});

dim = size(sig1);

sig1 = reshape(sig1, [], dim(3));
sig2 = reshape(sig2, [], dim(3));

if isempty(brain_mask)
    nanIdx = true(dim(1) * dim(2), 1);
else
    nanIdx = ~isnan(brain_mask(:));
end

sig1 = sig1(nanIdx, :)';
sig2 = sig2(nanIdx, :)';

if p.Results.detrend
    sig1 = detrend(sig1, 1);
    sig2 = detrend(sig2, 1);
end

r = NaN(maxlag * 2 + 1, dim(1) * dim(2));
r(:, nanIdx) = f_xcorr(sig1, sig2, maxlag);

r = reshape(r', dim(1), dim(2), maxlag * 2 + 1);

lag = -maxlag : maxlag;