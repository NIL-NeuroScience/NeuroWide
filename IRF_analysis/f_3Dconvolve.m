%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             f_3Dconvolve
% author - Brad Rauscher (created 2024)
% 
% Convolves each pixel of 'sig' by 'kernel' using the timing parameters in 
% 'win'. Only returns the times aligned to each frame in 'sig'.
% 
% INPUTS: f_3Dconvolve(sig, kernel, win, brain_mask)
%   sig: signal to convolve of size H x W x T
%   kernel: convolution kernel
%   win: timing parameters [t1, t2]
%   brain_mask: NaN mask of pixels to convolve of size H x W
% 
% OUTPUTS:
%   A: convolution product
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function A = f_3Dconvolve(sig, kernel, win, brain_mask)

nanIdx = isnan(brain_mask(:));

% reshape 'sig' to be T x numel(nanIdx)
dim = size(sig);
sig = reshape(sig, [], dim(3));
sig = sig(~nanIdx, :)';

% convolve 'sig' with 'kernel' and remove excess
A = conv2(sig, kernel);
A = A(-win(1) + 1 : end - win(2), :);

% reshape A back into size(sig)
tmp = nan(dim(1) * dim(2), dim(3));
tmp(~nanIdx, :) = A';

A = reshape(tmp, dim);