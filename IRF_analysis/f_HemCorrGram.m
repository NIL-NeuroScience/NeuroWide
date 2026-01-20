%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            f_HemCorrGram
% author - Brad Rauscher (created 2024)
% 
% Calculates correlation (pearson's coefficient) between 'sig1' and 'sig2'
% using a sliding window with parameters defined in 'win'.
% 
% INPUTS: f_HemCorrGram(sig1, sig2, win)
%   sig1: first signal
%   sig2: second signal
%   win: sliding time window [window size, slide increment]
% 
% OUTPUTS:
%   r: correlation value(s)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function r = f_HemCorrGram(sig1, sig2, win)

dim = size(sig1);

idx = (1 : win(2) : dim(3) - win(1) + 1)' + (0 : win(1) - 1);

frames = size(idx, 1);
tSig1 = zeros(dim(1) * frames, dim(2), win(1));
tSig2 = zeros(dim(1) * frames, dim(2), win(1));

for i = 1 : frames
    tSig1((i - 1) * dim(1) + 1 : i * dim(1), :, :) = sig1(:, :, idx(i, :));
    tSig2((i - 1) * dim(1) + 1 : i * dim(1), :, :) = sig2(:, :, idx(i, :));
end

std_tSig1 = std(tSig1, 0, 3);
std_tSig2 = std(tSig2, 0, 3);

tSig1 = tSig1 - mean(tSig1, 3);
tSig2 = tSig2 - mean(tSig2, 3);

calc_cov = (1 / win(1)) * sum(tSig1 .* tSig2, 3);

r = calc_cov ./ (std_tSig1 .* std_tSig2);
r = r';
r = reshape(r, [dim(2), dim(1), frames]);
r = permute(r, [2, 1, 3]);

end