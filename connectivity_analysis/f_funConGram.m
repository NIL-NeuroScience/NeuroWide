%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           f_funConGram
% author - Brad Rauscher (created 2024)
% 
% Calculates functional connectivity between columns of 'sig' over time
% window 'win'.
% 
% INPUTS: f_funConGram(sig, win)
%   sig: data in the form T (time) x N (channels)
%   win: sliding time window [window size, slide increment]
% 
% OUTPUTS:
%   r: correlation value(s)
%   f: figure handle
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [r, t] = f_funConGram(sig, win)

dim = size(sig);

% calculate starting index values for each window
% idx = 1 : win(2) : dim(1);
% idx(idx - 1 + win(1) > dim(1)) = [];
idx = (1 : win(2) : dim(1) - win(1) + 1)' + (0 : win(1) - 1);

% calculate functional connectivity of sig for each window
r = zeros(size(sig, 2), size(sig, 2), size(idx, 1));

for i = 1 : size(idx, 1)
    r(:, :, i) = corrcoef(sig(idx(i, :), :));
end

% calculate time values for each idx
t = win(1) / 2 : win(2) : dim(1) - win(1) / 2;