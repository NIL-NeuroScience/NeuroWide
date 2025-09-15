function [C,phi,fr] = f_hemCoherence(sig1,sig2,fs,mask,tapers)

%%
% sig1 = data.HbT;
% sig2 = data.gfp_HD;
% fs = 10;
% tapers = [5 9];
%

dim = size(sig1);
nanIdx = isnan(mask);
sig1 = reshape(sig1,dim(1)*dim(2),[]);
sig2 = reshape(sig2,dim(1)*dim(2),[]);
sig1 = sig1(~nanIdx,:)';
sig2 = sig2(~nanIdx,:)';

params = struct;
params.Fs = fs;
params.tapers = tapers;
params.trialave = 0;

[C,phi,~,~,~,fr] = coherencyc(sig1,sig2,params);

C = mean(C,2);
phi = mean(phi,2);
fr = fr';