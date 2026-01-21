%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           f_LR_varWeights
% author - Brad Rauscher (created 2024)
% 
% Performs optimization of linear regression weights and lags to predict 
% HbT from Ca and NE. Allows for spatially variant weights (A, B) and
% spatially invariant lags (t1, t2).
% 
% INPUTS: f_LR_varWeights(HbT, Ca, NE, win, brain_mask, _)
%   HbT: HbT matrix
%   Ca: Ca matrix
%   NE: NE matrix
%   win: window to use for IRF kernel (frames) [t1 t2]
%   brain_mask: mask of brain exposure (2D NaN image). Leave empty if no
%       mask is needed
% 
% OPTIONAL INPUTS:
%   ds: downsampling kernel (int; default 1)
%   maxThreads: maximum number of cores to use (default 4)
%   initialParams: initial timing parameters to use (default defined 
%       below)
%   corrWin: window for performance over time [t1 t2]. Leave empty to not
%       perform this analysis (default)
%   method: deconvolution method to use: 'direct' or 'fft' (default)
%   filter: band-pass filters Ca and NE (default = [0, 0.5])
%   norm: normalizes HbT, Ca, and NE matrices
% 
% OUTPUTS:
%   r: correlation map showing correlation between input HbT and
%       predicted HbT
%   outParams: output parameters including t1, t2, A, and B
%   conv_Ca: convolution product for Ca
%   conv_NE: convolution product for NE
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [r, outParams, conv_Ca, conv_NE] = ...
    f_LR_varWeights(HbT, Ca, NE, win, brain_mask, varargin)

p = inputParser;
addParameter(p, 'ds', 1);
addParameter(p, 'corrWin', []);
addParameter(p, 'maxThreads', 4);
addParameter(p, 'initialParams', 10 * [0.9, 0.1]);
addParameter(p, 'method', 'fft');
addParameter(p, 'filter', [0, 0.5]);
addParameter(p, 'norm', true);

parse(p, varargin{:});

maxNumCompThreads(p.Results.maxThreads);

if isempty(brain_mask)
    brain_mask = ones(size(Ca, [1, 2]));
end

if ~ismember(p.Results.method, {'fft', 'direct'})
    error("'method' value must be 'direct' or 'fft'.");
end

t1 = p.Results.initialParams(1);
t2 = p.Results.initialParams(2);

% downsample and reshape data

if ~isempty(p.Results.filter)
    Ca = f_bpf(Ca, p.Results.filter, 10, 3);
    NE = f_bpf(NE, p.Results.filter, 10, 3);
end

dsHbT = f_downsample(HbT, p.Results.ds);
dsCa = f_downsample(Ca, p.Results.ds);
dsNE = f_downsample(NE, p.Results.ds);
ds_brain_mask = f_downsample(brain_mask, p.Results.ds);

nanIdx = ~isnan(ds_brain_mask(:));
N = sum(nanIdx);

T = size(dsHbT, 3);
dsHbT = reshape(dsHbT, [], T);
dsCa = reshape(dsCa, [], T);
dsNE = reshape(dsNE, [], T);

dsHbT = dsHbT(nanIdx, :)';
dsCa = dsCa(nanIdx, :)';
dsNE = dsNE(nanIdx, :)';

% normalize data
if p.Results.norm
    dsHbT = dsHbT ./ std(dsHbT, 0, 1);
    dsCa = dsCa ./ std(dsCa, 0, 1);
    dsNE = dsNE ./ std(dsNE, 0, 1);
end

if string(p.Results.method) == "direct"
    % create design matrices
    l_irf = range(win) + 1;

    Ca_mat = f_2Dconvmtx(dsCa, l_irf);
    Ca_mat = Ca_mat(1 - win(1) : end - win(2), :, :);
    Ca_mat = permute(Ca_mat, [2, 3, 1]); % n_IRF x N x T

    NE_mat = f_2Dconvmtx(dsNE, l_irf);
    NE_mat = NE_mat(1 - win(1) : end - win(2), :, :);
    NE_mat = permute(NE_mat, [2, 3, 1]); % n_IRF x N x T
    
    design_matrix = [Ca_mat; NE_mat];
    design_HbT = permute(dsHbT, [3, 2, 1]);
    
    % optimize initial parameters
    
    A = ones(1, N);
    B = ones(1, N);
    
    hrf1 = f_LR_IRF(t1, A, win);
    hrf2 = f_LR_IRF(t2, B, win);
    
    conv_Ca = sum(Ca_mat .* hrf1, 1);
    conv_NE = sum(NE_mat .* hrf2, 1);
    
    LR = f_hemRegress(design_HbT, cat(4, conv_Ca, conv_NE), ones(1, N));
    
    A = LR(:, :, 1);
    B = LR(:, :, 2);
    
    % run optimization using gradient descent
    
    fun = @(params)f_hrf_cost_func(params(1), params(2), ...
        params(3 : N + 2), params(N + 3 : 2 * N + 2), win, ...
        design_matrix, design_HbT);
    
    options = optimset(...
        MaxFunEvals = 25000, ...
        MaxIter = 500, ...
        Display = 'off', ...
        Algorithm = 'active-set');

    params = fmincon(fun, [t1, t2, A, B], [], [], [], [], ...
        [0, win(1), -Inf * ones(1, numel(A) * 2)], ...
        [win(2), win(2), Inf * ones(1, numel(A) * 2)], [], options);
else
    A = ones(1, N);
    B = -ones(1, N);
    
    hrf1 = f_LR_IRF(t1, A, win);
    hrf2 = f_LR_IRF(t2, B, win);
    
    fft_length = 2^nextpow2(size(dsCa, 1) + size(hrf1, 1));
    fft_signal1 = fft(dsCa, fft_length);
    fft_signal2 = fft(dsNE, fft_length);

    convCa = ifft(fft_signal1 .* fft(hrf1, fft_length), 'symmetric');
    convNE = ifft(fft_signal2 .* fft(hrf2, fft_length), 'symmetric');

    convCa = convCa(1 - win(1) : T - win(1), :);
    convNE = convNE(1 - win(1) : T - win(1), :);

    LR = f_hemRegress(permute(dsHbT, [3, 2, 1]), ...
        permute(cat(4, convCa, convNE), [3, 2, 1, 4]), ones(1, N));
    
    A = LR(:, :, 1);
    B = -LR(:, :, 2);
    
    fun = @(params)f_hrf_cost_func_fft(params(1), params(2), ...
        params(3 : N + 2), params(N + 3 : 2 * N + 2), win, fft_signal1, ...
        fft_signal2, dsHbT, fft_length);
    
    options = optimset( ...
        MaxFunEvals = 25000, ...
        MaxIter = 500, ...
        Display = 'off', ...
        Algorithm = 'active-set');
    params = fmincon(fun, [t1, t2, A, B], [], [], [], [], ...
        [0, win(1), -Inf * ones(1, numel(A) * 2)], ...
        [win(2), win(2), Inf * ones(1, numel(A) * 2)], [], options);
end

% pixel-wise LR

hrf1 = f_LR_IRF(params(1), 1, win);
hrf2 = f_LR_IRF(params(2), 1, win);

if p.Results.norm
    conv_Ca = f_3Dconvolve(Ca ./ std(Ca, 0, 3), hrf1, win, brain_mask);
    conv_NE = f_3Dconvolve(NE ./ std(NE, 0, 3), hrf2, win, brain_mask);
    LR = f_hemRegress(HbT ./ std(HbT, 0, 3), ...
        cat(4, conv_Ca, conv_NE), brain_mask);
else
    conv_Ca = f_3Dconvolve(Ca, hrf1, win, brain_mask);
    conv_NE = f_3Dconvolve(NE, hrf2, win, brain_mask);
    LR = f_hemRegress(HbT, cat(4, conv_Ca, conv_NE), brain_mask);
end

outParams = struct;
outParams.tA = params(1);
outParams.tB = params(2);
outParams.A = LR(:, :, 1);
outParams.B = LR(:, :, 2);

% calculate performance

conv_Ca = conv_Ca .* LR(:, :, 1);
conv_NE = conv_NE .* LR(:, :, 2);

pred_HbT = conv_Ca + conv_NE;

if p.Results.norm
    r = f_corr(HbT ./ std(HbT, 0, 3), pred_HbT, 3);
else
    r = f_corr(HbT, pred_HbT, 3);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EXTRA FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% cost function for IRF optimization using direct method (MSE)
function J = f_hrf_cost_func(t1, t2, A, B, range, X_mat, y)
    HRF1 = f_LR_IRF(t1, A, range);
    HRF2 = f_LR_IRF(t2, B, range);
    
    hrf = [HRF1; HRF2];
    conv_result = sum(X_mat .* hrf, 1);
    J = mean((y - conv_result).^2, 'all');
end

% cost function for IRF optimization using fft method (MSE)
function J = f_hrf_cost_func_fft(t1, t2, A, B, range, CaMat, NEMat, y, n)
    HRF1 = f_LR_IRF(t1, A, range);
    HRF2 = f_LR_IRF(t2, B, range);

    fft_kernel1 = fft(HRF1, n);
    fft_kernel2 = fft(HRF2, n);

    conv_result1 = ifft(CaMat .* fft_kernel1, 'symmetric');
    conv_result1 = conv_result1(1 - range(1) : size(y, 1) - range(1), :);
    conv_result2 = ifft(NEMat .* fft_kernel2, 'symmetric');
    conv_result2 = conv_result2(1 - range(1) : size(y, 1) - range(1), :);

    J = mean((y - conv_result1 - conv_result2).^2, 'all');
end

% delta function IRF definition
function IRF = f_LR_IRF(t, A, range)
    tr = (range(1) : range(2))';
    IRF = zeros(numel(tr), 1);
    
    if rem(t, 1) == 0
        IRF(tr == t) = 1;
    else
        err = rem(1 + rem(t, 1), 1);
        IRF(tr == ceil(t)) = err;
        IRF(tr == floor(t)) = 1 - err;
    end
    IRF = IRF .* A;
end

end