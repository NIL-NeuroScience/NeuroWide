function [spectra,fr] = f_spectra(sig,fs,tapers)

%%
% sig = data.rfp_HD;
% fs = 10;
% fpass = [0 5];
% tapers = [5 9];
%

fpass = [0,0.5*fs];

params.Fs = fs;
params.fpass = fpass;
params.trialave = 0;
params.tapers = tapers;

[spectra,fr] = mtspectrumc(sig,params);

fr = fr';