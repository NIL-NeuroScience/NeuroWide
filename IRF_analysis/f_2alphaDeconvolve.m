function [perf,IRF,outParams,conv_Ca,conv_NE] = f_2alphaDeconvolve(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Performs deconvolution of HbT by applying a single alpha function 
% convolution kernel to Ca and NE, individually. Allows for spatially
% variant weights (A, B) for both alpha functions and spatially invariant
% timing parameters (tA,tB,tauA,tauB).
% INPUTS:
%   HbT - HbT video
%   Ca - Ca video
%   NE - NE video
%   win - window to use for IRF kernel (s) [t1 t2]
%   fs - sampling frequency (Hz)
%   brain_mask - mask of brain exposure (2D NaN image)
%   ds - downsampling kernel (int, typically 32)
%   maxThreads - maximum number of cores to use (default 4)
%   initialParams - initial timing parameters to use (default defined 
%       below)
% OUTPUTS:
%   perf - correlation map showing correlation between input HbT and
%       predicted HbT
%   IRF - both estimated IRFs
%   outParams - output parameters including tA, tB, tauA, tauB and brain
%       maps of A and B
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
addParameter(p,'initialParam',10*[0.1,-1,0.5,0.53]);
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
tau1 = p.Results.initialParam(3);
tau2 = p.Results.initialParam(4);

%% downsample and reshape data
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
    
    hrf1 = f_alpha_IRF(t1,tau1,0.5,A,0,win);
    hrf2 = f_alpha_IRF(t2,tau2,0.5,B,0,win);
    
    conv_Ca = sum(Ca_mat.*hrf1,1);
    conv_NE = sum(NE_mat.*hrf2,1);
    
    LR = f_hemRegress(design_HbT,cat(4,conv_Ca,conv_NE),ones(1,N));
    
    A = LR(:,:,1);
    B = LR(:,:,2);
    
    % run optimization using gradient descent
    
    fun = @(params)f_hrf_cost_func(params(1),params(2),params(3),params(4),params(5:N+4),params(N+5:2*N+4),win,design_matrix,design_HbT);
    
    options = optimset('MaxFunEvals',25000,'MaxIter',500,'Display','iter','Algorithm','active-set','Display','off');
    params = fmincon(fun,[t1, t2, tau1, tau2, A, B],[],[],[],[],[0 win(1) 0.01 0.01 -Inf*ones(1,numel(A)*2)],[win(2) win(2) win(2) win(2) Inf*ones(1,numel(A)*2)],[],options);
else
    A = ones(1,N);
    B = -ones(1,N);
    
    hrf1 = f_alpha_IRF(t1,tau1,0.5,A,0,win);
    hrf2 = f_alpha_IRF(t2,tau2,0.5,B,0,win);
    
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

    fun = @(params)f_hrf_cost_func_fft(params(1),params(2),params(3),params(4),params(5:N+4),params(N+5:2*N+4),win,fft_signal1,fft_signal2,dsHbT,fft_length);
    
    options = optimset('MaxFunEvals',25000,'MaxIter',500,'Display','iter','Algorithm','active-set','Display','off');
    params = fmincon(fun,[t1, t2, tau1, tau2, A, B],[],[],[],[],[0 win(1) 0.01 0.01 -Inf*ones(1,numel(A)*2)],[win(2) win(2) win(2) win(2) Inf*ones(1,numel(A)*2)],[],options);
end

%% pixel-wise LR

hrf1 = f_alpha_IRF(params(1),params(3),0.5,1,0,win);
hrf2 = f_alpha_IRF(params(2),params(4),0.5,1,0,win);

conv_Ca = f_3Dconvolve(Ca./std(Ca,0,3),hrf1,win,brain_mask);
conv_NE = f_3Dconvolve(NE./std(NE,0,3),hrf2,win,brain_mask);

LR = f_hemRegress(HbT./std(HbT,0,3),cat(4,conv_Ca,conv_NE),brain_mask);

outParams = struct;
outParams.tA = params(1);
outParams.tB = params(2);
outParams.tauA = params(3);
outParams.tauB = params(4);
outParams.A = LR(:,:,1);
outParams.B = LR(:,:,2);

IRF = [hrf1,-hrf2];

%% calculate performance

conv_Ca = conv_Ca.*LR(:,:,1);
conv_NE = conv_NE.*LR(:,:,2);

pred_HbT = conv_Ca+conv_NE;

perf = f_corr(HbT./std(HbT,0,3),pred_HbT,3);

%%

function J = f_hrf_cost_func(t1,t2,tau1,tau2,A,B,range,X_mat,y)
    HRF1 = f_alpha_IRF(t1,tau1,0.5,A,0,range);
    HRF2 = f_alpha_IRF(t2,tau2,0.5,B,0,range);

    hrf = [HRF1; HRF2];
    conv_result = sum(X_mat.*hrf,1);
    J = mean((y - conv_result).^2,'all');
end

function J = f_hrf_cost_func_fft(t1,t2,tau1,tau2,A,B,range,CaMat,NEMat,y,n)
    HRF1 = f_alpha_IRF(t1,tau1,0.5,A,0,range);
    HRF2 = f_alpha_IRF(t2,tau2,0.5,B,0,range);

    fft_kernel1 = fft(HRF1,n);
    fft_kernel2 = fft(HRF2,n);

    conv_result1 = ifft(CaMat.*fft_kernel1,'symmetric');
    conv_result1 = conv_result1(1-range(1):size(y,1)-range(1),:);
    conv_result2 = ifft(NEMat.*fft_kernel2,'symmetric');
    conv_result2 = conv_result2(1-range(1):size(y,1)-range(1),:);

    J = mean((y - conv_result1 - conv_result2).^2,'all');
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