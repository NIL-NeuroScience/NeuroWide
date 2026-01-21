%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            f_hemRegress
% author - Brad Rauscher (created 2024)
% 
% Solves the regression problem Y = A * X for each pixel in 'sig1'.
% 
% INPUTS: f_hemRegress(sig, reg, brain_mask)
%   Y: first signal
%   X: second signal. Can contain multiple regressors in dimension 4
%   brain_mask: mask of brain exposure (2D NaN image). Leave empty if no
%       mask is needed
% 
% OUTPUTS:
%   A: regression coefficients
%   residual: Y - A * X
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [A, residual] = f_hemRegress(Y, X, brain_mask)

% handle inputs
dim = size(Y);

if isempty(brain_mask)
    nanIdx = true(dim(1) * dim(2), 1);
else
    nanIdx = ~isnan(brain_mask(:));
end

% reshape X and Y
N_channels = size(X, 4);

Y = reshape(Y, dim(1) * dim(2), dim(3));
X = reshape(X, dim(1) * dim(2), dim(3), N_channels);

Y = Y(nanIdx, :)';
X = X(nanIdx, :, :);
X = permute(X, [2, 3, 1]);

% estimate A
A = zeros(sum(nanIdx), N_channels);

for i = 1 : sum(nanIdx)
    A(i, :) = X(:, :, i) \ Y(:, i);
end

% calculate residual and reshape outputs to original dimensions
residual = squeeze(sum(X .* permute(A, [3, 2, 1]), 2));
residual = Y - residual;

filler = NaN(dim(1) * dim(2), dim(3));
filler(nanIdx, :) = residual';

residual = reshape(filler, dim(1), dim(2), dim(3));

filler = NaN(dim(1) * dim(2), N_channels);
filler(nanIdx, :) = A;
A = reshape(filler, dim(1), dim(2), N_channels);