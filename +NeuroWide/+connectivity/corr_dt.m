%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               corr_dt
% author - Brad Rauscher (created 2023)
% 
% Calculates sliding correlation (pearson's coefficient) between matrices 
% sig1 and sig2 along 'dim' dimension. Removes average before calculating 
% correlation. Option to plot result as a correlation map.
% 
% INPUTS: f_corr(sig1, sig2, dim, _)
%   sig1: first signal
%   sig2: second signal
%   win: sliding window [width, shift]
%   dim: dimension to calculate correlation on
% 
% OUTPUTS:
%   r: correlation value(s)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function r = corr_dt(sig1, sig2, win, dim)

% check size
if size(sig1, dim) == size(sig2, dim)
    N = size(sig1, dim);
else
    error('Size of sig1 and sig2 do not agree!');
end

idx = 1:win(2):N - win(1) + 1;

idx = idx' + (0:win(1) - 1);

sig1_idx = repmat({':'}, 1, ndims(sig1));
sig2_idx = repmat({':'}, 1, ndims(sig2));

W = size(idx, 1);

r = cell(W, 1);

for i = 1:W
    sig1_idx{dim} = idx(i,:);
    sig2_idx{dim} = idx(i,:);
    
    reslice_1 = sig1(sig1_idx{:});
    reslice_2 = sig2(sig2_idx{:});
    
    % remove mean from sig1 and sig2
    reslice_1 = reslice_1 - mean(reslice_1, dim);
    reslice_2 = reslice_2 - mean(reslice_2, dim);
    
    % calculate pearson's coefficient along 'dim'
    std1 = sum(reslice_1.^2, dim);
    std2 = sum(reslice_2.^2, dim);
    std = sqrt(std1 .* std2);
    cov = sum((reslice_1 .* reslice_2), dim);
    
    r{i} = cov ./ std;
end

r = cat(dim, r{:});

end