function [perf,outParams,conv_Ca,conv_NE] = f_LR_varWeights(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Performs optimization of linear regression weights and lags to predict 
% HbT from Ca and NE. Allows for spatially variant weights (A, B) and
% spatially invariant lags (t1,t2).
% 
% INPUTS: f_1xIRF_varWeights(HbT,Ca,NE,win,brain_mask,_)
%   HbT - HbT video
%   Ca - Ca video
%   NE - NE video
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
%   outParams - output parameters including t0, tauA, tauB, A, and B
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                                                                                                                    
HbT = varargin{1};
Ca = varargin{2};
NE = varargin{3};
win = varargin{4};
brain_mask = varargin{5};

p = inputParser;
addParameter(p,'ds',1);
addParameter(p,'corrWin',[]);
addParameter(p,'maxThreads',4);
addParameter(p,'initialParam',10*[0.9,0.1]);
addParameter(p,'method','fft');

parse(p,varargin{6:end});

maxNumCompThreads(p.Results.maxThreads);

if isempty(brain_mask)
    brain_mask = ones(size(Ca,[1,2]));
end

if ~ismember(p.Results.method,{'fft','direct'})
    error("'method' value must be 'direct' or 'fft'.");
end

t1 = p.Results.initialParam(1);
t2 = p.Results.initialParam(2);

%% downsample and reshape data

Ca = f_bpf(Ca,[0,0.5],10,3);
NE = f_bpf(NE,[0,0.5],10,3);

dsHbT = f_downsample(HbT,p.Results.ds);
dsCa = f_downsample(Ca,p.Results.ds);
dsNE = f_downsample(NE,p.Results.ds);
ds_brain_mask = f_downsample(brain_mask,p.Results.ds);

nanIdx = ~isnan(ds_brain_mask(:));
N = sum(nanIdx);

[H,W,T] = size(dsHbT);
dsHbT = reshape(dsHbT,[],T);
dsCa = reshape(dsCa,[],T);
dsNE = reshape(dsNE,[],T);

dsHbT = dsHbT(nanIdx,:)';
dsCa = dsCa(nanIdx,:)';
dsNE = dsNE(nanIdx,:)';

% noramlize data
dsHbT = dsHbT./std(dsHbT,0,1);
dsCa = dsCa./std(dsCa,0,1);
dsNE = dsNE./std(dsNE,0,1);

if string(p.Results.method) == "direct"
    % create design matrices
    T = size(dsCa,1);
    l_irf = range(win)+1;
    idx_irf = win(1):win(2);
    
    i1 = abs(min([idx_irf;zeros(1,numel(idx_irf))],[],1))+1;
    i2 = T-max([idx_irf;zeros(1,numel(idx_irf))],[],1);
    i3 = max([idx_irf;zeros(1,numel(idx_irf))],[],1)+1;
    i4 = T-i1+1;
    
    Ca_mat = zeros(T,N,l_irf);
    NE_mat = zeros(T,N,l_irf);
    
    for v = 1:l_irf
        Ca_mat(i3(v):i4(v),:,v) = dsCa(i1(v):i2(v),:);
    end
    for v = 1:l_irf
        NE_mat(i3(v):i4(v),:,v) = dsNE(i1(v):i2(v),:);
    end
    
    Ca_mat = permute(Ca_mat,[3,2,1]);
    NE_mat = permute(NE_mat,[3,2,1]);
    
    design_matrix = [Ca_mat;NE_mat];
    design_HbT = permute(dsHbT,[3 2 1]);
    
    % optimize initial parameters
    
    A = ones(1,N);
    B = ones(1,N);
    
    hrf1 = f_LR_IRF(t1,A,win);
    hrf2 = f_LR_IRF(t2,B,win);
    
    conv_Ca = sum(Ca_mat.*hrf1,1);
    conv_NE = sum(NE_mat.*hrf2,1);
    
    LR = f_hemRegress(design_HbT,cat(4,conv_Ca,conv_NE),ones(1,N));
    
    A = LR(:,:,1);
    B = LR(:,:,2);
    
    % run optimization using gradient descent
    
    fun = @(params)f_hrf_cost_func(params(1),params(2),params(3:N+2),params(N+3:2*N+2),win,design_matrix,design_HbT);
    
    options = optimset('MaxFunEvals',25000,'MaxIter',500,'Display','iter','Algorithm','active-set','Display','off');
    params = fmincon(fun,[t1, t2, A, B],[],[],[],[],[0 win(1) -Inf*ones(1,numel(A)*2)],[win(2) win(2) Inf*ones(1,numel(A)*2)],[],options);
else
    A = ones(1,N);
    B = -ones(1,N);
    
    hrf1 = f_LR_IRF(t1,A,win);
    hrf2 = f_LR_IRF(t2,B,win);
    
    fft_length = 2^nextpow2(size(dsCa,1) + size(hrf1,1));
    fft_signal1 = fft(dsCa, fft_length);
    fft_signal2 = fft(dsNE, fft_length);

    convCa = ifft(fft_signal1 .* fft(hrf1,fft_length),'symmetric');
    convNE = ifft(fft_signal2 .* fft(hrf2,fft_length),'symmetric');

    convCa = convCa(1-win(1):T-win(1),:);
    convNE = convNE(1-win(1):T-win(1),:);

    LR = f_hemRegress(permute(dsHbT,[3,2,1]),permute(cat(4,convCa,convNE),[3,2,1,4]),ones(1,N));
    
    A = LR(:,:,1);
    B = -LR(:,:,2);
    
    fun = @(params)f_hrf_cost_func_fft(params(1),params(2),params(3:N+2),params(N+3:2*N+2),win,fft_signal1,fft_signal2,dsHbT,fft_length);
    
    options = optimset('MaxFunEvals',25000,'MaxIter',500,'Display','iter','Algorithm','active-set','Display','off');
    params = fmincon(fun,[t1, t2, A, B],[],[],[],[],[0 win(1) -Inf*ones(1,numel(A)*2)],[win(2) win(2) Inf*ones(1,numel(A)*2)],[],options);
end

%% pixel-wise LR

hrf1 = f_LR_IRF(params(1),1,win);
hrf2 = f_LR_IRF(params(2),1,win);

conv_Ca = f_3Dconvolve(Ca./std(Ca,0,3),hrf1,win,brain_mask);
conv_NE = f_3Dconvolve(NE./std(NE,0,3),hrf2,win,brain_mask);

LR = f_hemRegress(HbT./std(HbT,0,3),cat(4,conv_Ca,conv_NE),brain_mask);

outParams = struct;
outParams.tA = params(1);
outParams.tB = params(2);
outParams.A = LR(:,:,1);
outParams.B = LR(:,:,2);

%% calculate performance

conv_Ca = conv_Ca.*LR(:,:,1);
conv_NE = conv_NE.*LR(:,:,2);

pred_HbT = conv_Ca+conv_NE;

perf = f_corr(HbT./std(HbT,0,3),pred_HbT,3);

%%

function J = f_hrf_cost_func(t1,t2,A,B,range,X_mat,y)
    HRF1 = f_LR_IRF(t1,A,range);
    HRF2 = f_LR_IRF(t2,B,range);
    
    hrf = [HRF1; HRF2];
    conv_result = sum(X_mat.*hrf,1);
    J = mean((y - conv_result).^2,'all');
end

function J = f_hrf_cost_func_fft(t1,t2,A,B,range,CaMat,NEMat,y,n)
    HRF1 = f_LR_IRF(t1,A,range);
    HRF2 = f_LR_IRF(t2,B,range);

    fft_kernel1 = fft(HRF1,n);
    fft_kernel2 = fft(HRF2,n);

    conv_result1 = ifft(CaMat.*fft_kernel1,'symmetric');
    conv_result1 = conv_result1(1-range(1):size(y,1)-range(1),:);
    conv_result2 = ifft(NEMat.*fft_kernel2,'symmetric');
    conv_result2 = conv_result2(1-range(1):size(y,1)-range(1),:);

    J = mean((y - conv_result1 - conv_result2).^2,'all');
end

function IRF = f_LR_IRF(t,A,range)
    hrf_l = range;

    tr = (hrf_l(1):hrf_l(2))';
    IRF = zeros(numel(tr),1);
    
    if rem(t,1) == 0
        IRF(tr==t) = 1;
    else
        err = rem(1+rem(t,1),1);
        IRF(tr==ceil(t)) = err;
        IRF(tr==floor(t)) = 1-err;
    end
    IRF = IRF.*A;
end

end