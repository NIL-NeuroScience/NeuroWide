function [conv_sig] = f_3Dconvolve(sig,kernel,win,brain_mask)

nanIdx = isnan(brain_mask(:));

dim = size(sig);
sig = reshape(sig,[],dim(3));
sig = sig(~nanIdx,:)';

conv_sig = conv2(sig,kernel);
conv_sig = conv_sig(-win(1)+1:end-win(2),:);

tmp = nan(dim(1)*dim(2),dim(3));
tmp(~nanIdx,:) = conv_sig';

conv_sig = reshape(tmp,dim(1),dim(2),dim(3));