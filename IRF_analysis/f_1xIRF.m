function [perf,IRF,outParams,corrGram] = f_1xIRF(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Performs deconvolution of HbT by applying a single double alpha function 
% convolution kernel to Ca.
% 
% INPUTS: f_1xIRF(HbT,Ca,win,brain_mask,_)
%   HbT - HbT video
%   Ca - Ca video
%   win - window to use for IRF kernel (frames) [t1 t2]
%   brain_mask - mask of brain exposure (2D NaN image). Leave empty if no
%       mask is needed
% 
% EXTRA PARAMETERS:
%   ds - downsampling kernel (int; default 1)
%   maxThreads - maximum number of cores to use (default 4)
%   initialParams - initial timing parameters to use (default defined 
%       below)
%   corrWin - window for performance over time [t1 t2]. Leave empty to not
%       perform this analysis (default)
%   method - deconvolution method to use: 'direct' or 'fft' (default)
% 
% OUTPUTS:
%   perf - correlation map showing correlation between input HbT and
%       predicted HbT
%   IRF - estimated IRFs across the cortex
%   outParams - output parameters including t0, tauA, tauB, A, and B
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

HbT = varargin{1};
Ca = varargin{2};
win = varargin{3};
brain_mask = varargin{4};

p = inputParser;
addParameter(p,'ds',1);
addParameter(p,'corrWin',[]);
addParameter(p,'maxThreads',4);
addParameter(p,'initialParam',10*[0.1,0.5,0.53]);
addParameter(p,'method','fft');

parse(p,varargin{5:end});

maxNumCompThreads(p.Results.maxThreads);

if isempty(brain_mask)
    brain_mask = ones(size(Ca,[1,2]));
end

if ~ismember(p.Results.method,{'fft','direct'})
    error("'method' value must be 'direct' or 'fft'.");
end

t0 = p.Results.initialParam(1);
tau1 = p.Results.initialParam(2);
tau2 = p.Results.initialParam(3);

%% downsample and reshape data
dsHbT = f_downsample(HbT,p.Results.ds);
dsCa = f_downsample(Ca,p.Results.ds);
ds_brain_mask = f_downsample(brain_mask,p.Results.ds);

nanIdx = ~isnan(ds_brain_mask(:));
N = sum(nanIdx);

[H,W,T] = size(dsHbT);
dsHbT = reshape(dsHbT,[],T);
dsCa = reshape(dsCa,[],T);

dsHbT = dsHbT(nanIdx,:)';
dsCa = dsCa(nanIdx,:)';

%% normalize data
dsHbT = dsHbT./std(dsHbT,0,1);
dsCa = dsCa./std(dsCa,0,1);

%% optimization algorithm
if string(p.Results.method) == "direct"
    % create design matrices
    
    n_IRF = (win(2)-win(1))+1;
    
    Ca_mat = zeros(T,N,n_IRF);
    
    T = size(dsCa,1);
    l_irf = range(win)+1;
    idx_irf = win(1):win(2);
    
    i1 = abs(min([idx_irf;zeros(1,numel(idx_irf))],[],1))+1;
    i2 = T-max([idx_irf;zeros(1,numel(idx_irf))],[],1);
    i3 = max([idx_irf;zeros(1,numel(idx_irf))],[],1)+1;
    i4 = T-i1+1;
    
    for v = 1:l_irf
        Ca_mat(i3(v):i4(v),:,v) = dsCa(i1(v):i2(v),:);
    end
    
    Ca_mat = reshape(Ca_mat,[],n_IRF);
    design_HbT = dsHbT(:);
    
    % optimize initial parameters    
    A = 1;
    B = -1;
    
    hrf1 = f_alpha_IRF(t0,tau1,tau2,A,0,win);
    hrf2 = f_alpha_IRF(t0,tau1,tau2,0,B,win);
    
    convPos = sum(Ca_mat.*hrf1',2);
    convNeg = sum(Ca_mat.*hrf2',2);
    
    LR = [convPos,convNeg] \ design_HbT;
    A = LR(1);
    B = -LR(2);
    
    % run optimization using gradient descent
    fun = @(params)f_hrf_cost_func(params(1),params(2),params(3),params(4),params(5),win,Ca_mat,design_HbT);
    
    options = optimset('MaxFunEvals',25000,'MaxIter',500,'Display','iter','Algorithm','active-set','FunValCheck','on','Display','off');
    params = fmincon(fun,[t0, tau1, tau2, A, B],[],[],[],[],[0,0.01,0.01,-Inf,-Inf],[win(2),win(2),win(2),Inf,Inf],[],options);
else
    A = 1;
    B = -1;
    
    hrf1 = f_alpha_IRF(t0,tau1,tau2,A,0,win);
    hrf2 = f_alpha_IRF(t0,tau1,tau2,0,B,win);
    
    fft_length = 2^nextpow2(size(dsCa,1) + numel(hrf1));
    fft_signal = fft(dsCa, fft_length);

    convPos = ifft(fft_signal .* fft(hrf1,fft_length),'symmetric');
    convNeg = ifft(fft_signal .* fft(hrf2,fft_length),'symmetric');

    convPos = convPos(1:T,:);
    convNeg = convNeg(1:T,:);

    LR = [convPos(:),convNeg(:)] \ dsHbT(:);
    A = LR(1);
    B = -LR(2);
    
    %
    fun = @(params)f_hrf_cost_func_fft(params(1),params(2),params(3),params(4),params(5),win,fft_signal,dsHbT,fft_length);
    
    options = optimset('MaxFunEvals',25000,'MaxIter',500,'Display','iter','Algorithm','active-set','FunValCheck','on','Display','off');
    params = fmincon(fun,[t0, tau1, tau2, A, B],[],[],[],[],[0,0.01,0.01,-Inf,-Inf],[win(2),win(2),win(2),Inf,Inf],[],options);
end

%% calculate sliding correlation window between predicted HbT and HbT
IRF = f_alpha_IRF(params(1),params(2),params(3),params(4),params(5),win);
conv = f_3Dconvolve(Ca./std(Ca,0,3),IRF,win,ones(size(brain_mask)));
perf = f_corr(HbT,conv,3);

if ~isempty(p.Results.corrWin)
    corrGram = f_HemCorrGram(HbT,conv,p.Results.corrWin);
end

%% save params
outParams = struct;
outParams.t0 = params(1);
outParams.tau1 = params(2);
outParams.tau2 = params(3);
outParams.A = params(4);
outParams.B = params(5);

%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function J = f_hrf_cost_func(t0,tau1,tau2,A,B,range,X_mat,y)
    hrf = f_alpha_IRF(t0,tau1,tau2,A,B,range);
    
    conv_result = sum(X_mat.*hrf',2);
    J = mean((y - conv_result).^2,'all');
end

function J = f_hrf_cost_func_fft(t0,tau1,tau2,A,B,range,X_mat,y,n)
    hrf = f_alpha_IRF(t0,tau1,tau2,A,B,range);
    fft_kernel = fft(hrf,n);

    conv_result = ifft(X_mat.*fft_kernel,'symmetric');
    conv_result = conv_result(1:size(y,1),:);

    J = mean((y - conv_result).^2,'all');
end

function IRF = f_alpha_IRF(t0,tau1,tau2,A,B,range)
    hrf_l = range;
    tr = (((hrf_l(1):hrf_l(2)))-t0)';
    D = ((tr)./tau1).^3.*exp(-(tr)./tau1);
    D(tr<0) = 0;
    C = ((tr)./tau2).^3.*exp(-(tr)./tau2);
    C(tr<0) = 0;
    
    IRF = A.*D + B.*C;
end

end
