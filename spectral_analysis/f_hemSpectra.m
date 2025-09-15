function [spectra,fr] = f_hemSpectra(sig,fs,fpass,mask,tapers)

%%
% sig = data.rfp_HD;
% fs = 10;
% fpass = [0 5];
% tapers = [5 9];
%

dim = size(sig);
sig = reshape(sig,[dim(1)*dim(2), dim(3)]);
nanIdx = isnan(mask);
sig = sig(~nanIdx,:);
sig = sig';

params.Fs = fs;
params.fpass = fpass;
params.trialave = 0;
params.tapers = tapers;

[spectraMat,fr] = mtspectrumc(sig,params);

spectra = mean(spectraMat,2);
fr = fr';



