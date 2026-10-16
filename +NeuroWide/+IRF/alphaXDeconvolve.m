%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          alpha2Deconvolve
% author - Brad Rauscher (created 2024)
% 
% Performs deconvolution of HbT by applying a single alpha function 
% convolution kernel to Ca and NE, individually. Allows for spatially
% variant weights (A, B) for both alpha functions and spatially invariant
% timing parameters (tA,tB,tauA,tauB).
% 
% INPUTS: f_2alphaDeconvolve(Y, X, win, _)
%   Y: regressee matrix
%   X: regressor matrix
%   win: window to use for IRF kernel (frames) [t1 t2]
% 
% OPTIONAL INPUTS:
%   ds: downsampling kernel (int; default 1)
%   maxThreads: maximum number of cores to use (default 4)
%   initialParams: initial timing parameters to use (default defined 
%       below)
%   corrWin: window for performance over time [t1 t2]. Leave empty to not
%       perform this analysis (default)
%   method: deconvolution method to use: 'direct' or 'fft' (default)
% 
% OUTPUTS:
%   r: correlation map showing correlation between input HbT and
%       predicted HbT
%   IRF: estimated IRFs across the cortex
%   outParams: output parameters including t0, tauA, tauB, A, and B
%   conv_Ca: convolution product of IRF_Ca with Ca
%   conv_NE: convolution product of IRF_NE with NE
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [r, IRF, outParams, conv] = ...
    alphaXDeconvolve(Y, X, win, varargin)

%%
p = inputParser;
addParameter(p, 'corrWin', []);
addParameter(p, 'maxThreads', 4);
addParameter(p, 'initialParam', 10 * [0.1, -1, 0.5, 0.53]);

parse(p,varargin{:});

maxNumCompThreads(p.Results.maxThreads);

t1 = p.Results.initialParam(1);
t2 = p.Results.initialParam(2);
tau1 = p.Results.initialParam(3);
tau2 = p.Results.initialParam(4);

%% downsample and reshape data

[T, C, N] = size(X);

% noramlize data
normY = Y - mean(Y);
normX = X - mean(X);
normY = normY ./ std(normY, 0, 1);
normX = normX ./ std(normX, 0, 1);

%%
A = ones(1, C);

hrf = NeuroWide.IRF.alpha_IRF(t1, tau1, 0.5, A, 0, win);
l_irf = range(win) + 1;

fft_length = 2^nextpow2(T + l_irf);
fft_X = fft(normX, fft_length);

convX = ifft(fft_X .* fft(hrf, fft_length), 'symmetric');

convX = convX(1 - win(1) : T - win(1), :, :);

LR = multiRegress(normY, convX);

initial_params = [t1 * ones(N, 1);
                    tau1 * ones(N, 1);
                    LR(:)];
lower_bound = [win(1) * ones(N, 1);
                    0.01 * ones(N, 1);
                    -Inf * ones(numel(LR), 1)];
upper_bound = [win(2) * ones(N, 1);
                    win(2) * ones(N, 1);
                    Inf * ones(numel(LR), 1)];

%%
fun = @(params)f_hrf_cost_func_fft(params(1:N), params(N+1:2*N), ...
    params(2*N+1:end), win, ...
    fft_X, normY, fft_length);

options = optimset( ...
    MaxFunEvals = 25000, ...
    MaxIter = 500, ...
    Display = 'off', ...
    Algorithm = 'active-set');

params = fmincon(fun, initial_params, [], [], [], [], ...
    lower_bound, upper_bound, ...
    [], options);

%% pixel-wise LR

t = params(1:N);
tau = params(N+1:2*N);

hrf = zeros(range(win)+1, N);
for i = 1:N
    hrf(:,i) = NeuroWide.IRF.alpha_IRF(t(i), tau(i), 0.5, 1, 0, win);
end

conv = zeros(T, C, N);
for i = 1:N
    tmp = conv2(normX(:,:,i), hrf(:,i));
    conv(:,:,i) = tmp(-win(1) + 1 : end - win(2), :);
end

LR = multiRegress(normY, conv);

%%
outParams = struct;
outParams.t = t;
outParams.tau = tau;
outParams.A = LR;

IRF = hrf;

% calculate performance

conv = conv .* permute(LR,[3,1,2]);

r = NeuroWide.connectivity.corr(normY, sum(conv,3), 1);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EXTRA FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% cost function for IRF optimization using direct method (MSE)
function J = f_hrf_cost_func(t1, t2, tau1, tau2, A, B, range, X_mat, y)
    HRF1 = NeuroWide.IRF.alpha_IRF(t1, tau1, 0.5, A, 0, range);
    HRF2 = NeuroWide.IRF.alpha_IRF(t2, tau2, 0.5, B, 0, range);

    hrf = [HRF1; HRF2];
    conv_result = sum(X_mat .* hrf, 1);
    J = mean((y - conv_result).^2, 'all');
end

% cost function for IRF optimization using fft method (MSE)
function J = f_hrf_cost_func_fft(t, tau, A, win, ...
        data, y, fft_length)
    
    N = numel(t);
    A = reshape(A, [], N);
    C = size(A, 1);

    IRF = zeros(range(win)+1, C, N);
    for i = 1:N
        IRF(:,:,i) = NeuroWide.IRF.alpha_IRF(t(i), tau(i), 0.5, A(:,i)', 0, win);
    end

    fft_kernel = fft(IRF, fft_length, 1);

    conv_result = ifft(data .* fft_kernel, 'symmetric');
    conv_result = conv_result(1 - win(1) : size(y, 1) - win(1), :, :);

    J = mean((y - sum(conv_result, 3)).^2, 'all');
end

function LR = multiRegress(Y, X)

    C = size(Y, 2);
    N = size(X, 3);

    LR = zeros(C, N);

    for i = 1:C
        LR(i,:) = squeeze(X(:,i,:)) \ Y(:,i);
    end

end

end