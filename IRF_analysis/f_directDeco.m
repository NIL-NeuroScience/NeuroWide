function [perf,IRF] = f_directDeco(toPred,predictor,win,fs, ...
    brain_mask,ds)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Performs deconvolution of HbT by applying a single alpha function 
% convolution kernel to Ca and NE, individually. Allows for spatially
% variant weights (A, B) for both alpha functions and spatially invariant
% timing parameters (tA,tB,tauA,tauB).
% INPUTS:
%   HbT - HbT video
%   Ca - Ca video
%   win - window to use for IRF kernel (s) [t1 t2]
%   fs - sampling frequency (Hz)
%   brain_mask - mask of brain exposure (2D NaN image)
%   ds - downsampling kernel (int, typically 32)
% OUTPUTS:
%   perf - correlation map showing correlation between input HbT and
%       predicted HbT
%   IRF - both estimated IRFs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% downsample and reshape data
dstoPred = f_downsample(toPred,ds);
dspredictor = f_downsample(predictor,ds);

ds_brain_mask = f_downsample(brain_mask,ds);
nanIdx = ~isnan(ds_brain_mask(:));
N = sum(nanIdx);

dim = size(dstoPred);
dstoPred = reshape(dstoPred,[],dim(3));
dspredictor = reshape(dspredictor,[],dim(3));

dstoPred = dstoPred(nanIdx,:)';
dspredictor = dspredictor(nanIdx,:)';

dstoPred = dstoPred./std(dstoPred,0,1);
dspredictor = dspredictor./std(dspredictor,0,1);

%% create design matrices
T = size(dspredictor,1);
l_irf = fs*range(win)+1;
idx_irf = win(1)*fs:win(2)*fs;

i1 = abs(min([idx_irf;zeros(1,numel(idx_irf))],[],1))+1;
i2 = T-max([idx_irf;zeros(1,numel(idx_irf))],[],1);
i3 = max([idx_irf;zeros(1,numel(idx_irf))],[],1)+1;
i4 = T-i1+1;

predictor_mat = zeros(T,N,l_irf);

for v = 1:l_irf
    predictor_mat(i3(v):i4(v),:,v) = dspredictor(i1(v):i2(v),:);
end

design_matrix = reshape(predictor_mat,[],l_irf);
design_toPred = dstoPred(:);

%%

IRF = design_matrix \ design_toPred;

%% calculate performance

pred_sig = f_3Dconvolve(predictor,IRF,win*fs,ones(size(brain_mask)));
perf = f_corr(toPred,pred_sig,3);

end