%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          f_1xIRF_varWeights
% author - Brad Rauscher (created 2024)
% 
% Performs deconvolution of HbT by applying a single double alpha function 
% convolution kernel to Ca. Allows for spatially variant weights (A, B) for
% both alpha functions and spatially invariant timing parameters 
% (t0,tauA,tauB).
% 
% INPUTS: f_1xIRF_varWeights(HbT, Ca, win, brain_mask, _)
%   HbT: HbT matrix
%   Ca: Ca matrix
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
%   norm: normalizes Ca and HbT matrices
% 
% OUTPUTS:
%   r: correlation map showing correlation between input HbT and
%       predicted HbT
%   IRF: estimated IRFs across the cortex
%   outParams: output parameters including t0, tauA, tauB, A, and B
%   r_dt: correlation matrix between 'HbT' and 'pred_HbT'
%   pred_HbT: predicted HbT values using 'IRF'
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [r, IRF, outParams, r_dt, pred_HbT] = ...
    f_1xIRF_varWeights(HbT, Ca, win, brain_mask, varargin)

% parse additional inputs
p = inputParser;
addParameter(p, 'ds' ,1); % downsampling factor (square binning)
addParameter(p, 'corrWin', []); % [window size, increment] to calculate 
% sliding window correlation between HbT and pred_HbT
addParameter(p, 'maxThreads', 4); % max number of threads
addParameter(p, 'initialParam', 10 * [0.1, 0.5, 0.53]); % initial 
% parameters
addParameter(p, 'method', 'fft'); % model for fitting parameters
addParameter(p, 'norm' ,1); % normalize HbT and Ca

parse(p, varargin{:});

% set max number of threads for fitting algorithm
maxNumCompThreads(p.Results.maxThreads);

% handle inputs
if isempty(brain_mask)
    brain_mask = ones(size(Ca, [1, 2]));
end

if ~ismember(p.Results.method, {'fft', 'direct'})
    error("'method' value must be 'direct' or 'fft'.");
end

t0 = p.Results.initialParam(1);
tau1 = p.Results.initialParam(2);
tau2 = p.Results.initialParam(3);

% downsample and reshape data
dsHbT = f_downsample(HbT, p.Results.ds);
dsCa = f_downsample(Ca, p.Results.ds);
ds_brain_mask = f_downsample(brain_mask, p.Results.ds);

nanIdx = ~isnan(ds_brain_mask(:));
N = sum(nanIdx);

T = size(dsHbT, 3);
dsHbT = reshape(dsHbT, [], T);
dsCa = reshape(dsCa, [], T);

dsHbT = dsHbT(nanIdx, :)';
dsCa = dsCa(nanIdx, :)';

% normalize data
if p.Results.norm
    dsHbT = dsHbT ./ std(dsHbT, 0, 1);
    dsCa = dsCa ./ std(dsCa, 0, 1);
end

% optimization algorithm
if string(p.Results.method) == "direct"
    % create design matrices
    n_IRF = (win(2) - win(1)) + 1;
    
    Ca_mat = f_2Dconvmtx(dsCa, n_IRF);
    Ca_mat = Ca_mat(1 - win(1) : end - win(2), :, :);
    Ca_mat = permute(Ca_mat, [2, 3, 1]); % n_IRF x N x T
    
    design_HbT = permute(dsHbT, [3, 2, 1]);
    
    % optimize initial parameters
    
    A = ones(1, N);
    B = -ones(1, N);
    
    hrf1 = f_alpha_IRF(t0, tau1, tau2, A, 0, win);
    hrf2 = f_alpha_IRF(t0, tau1, tau2, 0, B, win);
    
    convPos = sum(Ca_mat .* hrf1, 1);
    convNeg = sum(Ca_mat .* hrf2, 1);
    
    LR = f_hemRegress(design_HbT, cat(4, convPos, convNeg), ones(1, N));
    
    A = LR(:, :, 1);
    B = -LR(:, :, 2);
    
    % run optimization using gradient descent
    fun = @(params)f_hrf_cost_func(params(1), params(2), params(3), ...
        params(4 : N + 3), params(N + 4 : 2 * N + 3), win, Ca_mat, ...
        design_HbT);
    
    options = optimset(...
        MaxFunEvals = 25000, ...
        MaxIter = 500, ...
        Display = 'off', ...
        Algorithm = 'active-set');

    params = fmincon(fun, [t0, tau1, tau2, A, B], [], [], [], [], ...
        [0, 0.01, 0.01, -Inf * ones(1, numel(A) * 2)], [win(2), ...
        win(2), win(2), Inf * ones(1, numel(A) * 2)], [], options);

else
    % FIR FFT deconvolution
    A = ones(1, N);
    B = -ones(1, N);
    
    hrf1 = f_alpha_IRF(t0, tau1, tau2, A, 0, win);
    hrf2 = f_alpha_IRF(t0, tau1, tau2, 0, B, win);
    
    fft_length = 2^nextpow2(size(dsCa, 1) + size(hrf1, 1));
    fft_signal = fft(dsCa, fft_length);

    convPos = ifft(fft_signal .* fft(hrf1, fft_length), 'symmetric');
    convNeg = ifft(fft_signal .* fft(hrf2, fft_length), 'symmetric');

    convPos = convPos(1 : T, :);
    convNeg = convNeg(1 : T, :);

    LR = f_hemRegress(permute(dsHbT, [3, 2, 1]), ...
        permute(cat(4, convPos, convNeg), [3, 2, 1, 4]), ones(1, N));
    
    A = LR(:, :, 1);
    B = -LR(:, :, 2);

    % optimization algorithm
    fun = @(params)f_hrf_cost_func_fft(params(1), params(2), params(3), ...
        params(4 : N + 3), params(N + 4 : 2 * N + 3), win, fft_signal, ...
        dsHbT, fft_length);
    
    options = optimset('MaxFunEvals', 25000, ...
        MaxIter = 500, ...
        Display = 'off', ...
        Algorithm = 'active-set');
    
    params = fmincon(fun, [t0, tau1, tau2, A, B], [], [], [], [], ...
        [0, 0.01, 0.01, -Inf * ones(1, numel(A) * 2)], [win(2), ...
        win(2), win(2), Inf * ones(1, numel(A) * 2)], [], options);
end

% pixel-wise LR
hrf1 = f_alpha_IRF(params(1), params(2), params(3), 1, 0, win);
hrf2 = f_alpha_IRF(params(1), params(2), params(3), 0, -1, win);

if p.Results.norm
    normCa = Ca ./ std(Ca, 0, 3);
    normHbT = HbT ./ std(HbT, 0, 3);
else
    normCa = Ca;
    normHbT = HbT;
end
convPos = f_3Dconvolve(normCa, hrf1, win, brain_mask);
convNeg = f_3Dconvolve(normCa, hrf2, win, brain_mask);

LR = f_hemRegress(normHbT, cat(4, convPos, convNeg), brain_mask);

% save params
outParams = struct;
outParams.t0 = params(1);
outParams.tauA = params(2);
outParams.tauB = params(3);
outParams.A = LR(:, :, 1);
outParams.B = -LR(:, :, 2);

% recalculate IRF
dim = size(Ca);

IRF = f_alpha_IRF(params(1), params(2), params(3), outParams.A(:)', ...
    outParams.B(:)', win)';
IRF = reshape(IRF, dim(1), dim(2), []);

% calculate performance
pred_HbT = convPos .* LR(:, :, 1) + convNeg .* LR(:, :, 2);

r = f_corr(HbT, pred_HbT, 3);

if ~isempty(p.Results.corrWin)
    r_dt = f_HemCorrGram(HbT, pred_HbT, p.Results.corrWin);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EXTRA FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% cost function for IRF optimization using direct method (MSE)
function J = f_hrf_cost_func(t0, tau1, tau2, A, B, range, X_mat, y)
    hrf = f_alpha_IRF(t0, tau1, tau2, A, B, range);
    
    conv_result = sum(X_mat .* hrf, 1);
    J = mean((y - conv_result).^2, 'all');
end

% cost function for IRF optimization using fft method (MSE)
function J = f_hrf_cost_func_fft(t0, tau1, tau2, A, B, range, X_mat, y, n)
    hrf = f_alpha_IRF(t0, tau1, tau2, A, B, range);
    fft_kernel = fft(hrf, n);

    conv_result = ifft(X_mat .* fft_kernel, 'symmetric');
    conv_result = conv_result(1 : size(y, 1), :);

    J = mean((y - conv_result).^2, 'all');
end

% double gamma function definition
function IRF = f_alpha_IRF(t0, tau1, tau2, A, B, range)
    hrf_l = range;
    tr = ((hrf_l(1) : hrf_l(2)) - t0)';
    D = (tr ./tau1).^3 .* exp(-tr ./ tau1);
    D(tr<0) = 0;
    C = (tr ./ tau2).^3 .* exp(-tr ./ tau2);
    C(tr < 0) = 0;
    
    IRF = A .* D + B .* C;
end

end