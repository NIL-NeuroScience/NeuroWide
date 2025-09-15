function [xcorr,lag] = f_hemLag_dT(sig1,sig2,fs,lagwin,mask)

%%
% sig1 = data.HbT;
% sig2 = data.rfp_HD;
% fs = 10;
% lagwin = [-5 5];
%

dim = size(sig1);

sig1 = reshape(sig1,dim(1)*dim(2),dim(3));
sig2 = reshape(sig2,dim(1)*dim(2),dim(3));
nanIdx = isnan(mask);
sig1 = sig1(~nanIdx,:)';
sig2 = sig2(~nanIdx,:)';
sig1 = detrend(sig1,1);
sig2 = detrend(sig2,1);

maxlag = fs*max(abs(lagwin));
lagRange = -maxlag:maxlag;lagRange = lagRange >= lagwin(1)*fs & lagRange <= lagwin(2)*fs;

xcorrMat = f_xcorr(sig1,sig2,maxlag);
badPixels = isnan(sum(xcorrMat,1));

xcorr = mean(xcorrMat(:,~badPixels),2);
lag = (-maxlag:maxlag)/fs;