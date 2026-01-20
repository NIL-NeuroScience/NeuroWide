%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            f_directDeco
% author - Brad Rauscher (created 2024)
% 
% Performs direct deconvolution to estimate convolution kernel 'IRF' in
% X * IRF = Y.
% 
% INPUTS: f_directDeco(Y, X, win, brain_mask, _)
%   Y: signal to predict
%   X: predictor signal
%   win: window to use for IRF kernel (frames) [t1 t2]
%   brain_mask: mask of brain exposure (2D NaN image).
% 
% OPTIONAL INPUTS:
%   ds: downsampling kernel (int; default 1)
%   norm: normalizes X and Y matrices
% 
% OUTPUTS:
%   r: correlation map showing correlation between input Y and predicted Y
%   IRF: estimated IRF
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [r, IRF] = f_directDeco(Y, X, win, brain_mask, varargin)

% parse inputs
p = inputParser;
addParameter(p, 'ds', 1); % downsampling factor (square binning)
addParameter(p, 'norm', 1); % normalize HbT and Ca

parse(p, varargin{:});
ds = p.Results.ds;

% downsample and reshape data
dstoPred = f_downsample(Y, ds);
dspredictor = f_downsample(X, ds);

ds_brain_mask = f_downsample(brain_mask, ds);
nanIdx = ~isnan(ds_brain_mask(:));

dim = size(dstoPred);
dstoPred = reshape(dstoPred, [], dim(3));
dspredictor = reshape(dspredictor, [], dim(3));

dstoPred = dstoPred(nanIdx, :)';
dspredictor = dspredictor(nanIdx, :)';

if p.Results.norm
    dstoPred = dstoPred ./ std(dstoPred, 0, 1);
    dspredictor = dspredictor ./ std(dspredictor, 0, 1);
end

% create design matrices
l_irf = range(win) + 1;

predictor_mat = f_2Dconvmtx(dspredictor, l_irf);
predictor_mat = predictor_mat(1 - win(1) : end - win(2), :, :);
predictor_mat = permute(predictor_mat, [1, 3, 2]);

design_matrix = reshape(predictor_mat, [], l_irf);
design_toPred = dstoPred(:);

% estimate IRF
IRF = design_matrix \ design_toPred;

% calculate performance
pred_sig = f_3Dconvolve(X, IRF, win, ones(size(brain_mask)));
r = f_corr(Y, pred_sig, 3);

end