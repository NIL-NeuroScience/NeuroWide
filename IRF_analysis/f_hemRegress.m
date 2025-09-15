function [rMat,residual] = f_hemRegress(sig,reg,brain_mask)

dim = size(sig);

sig = reshape(sig,[],dim(3));
reg = reshape(reg,dim(1)*dim(2),dim(3),[]);
nanIdx = isnan(brain_mask(:));

sig = sig(~nanIdx,:);
sig = sig';
reg = reg(~nanIdx,:,:);
reg = permute(reg,[2 3 1]);

N = size(sig,2);

rMat = zeros(N,size(reg,2));

for idx = 1:N
    rMat(idx,:) = reg(:,:,idx) \ sig(:,idx);
end

residual = squeeze(sum(reg.*permute(rMat,[3 2 1]),2));
residual = (sig-residual)';

nanMat = NaN(dim(1)*dim(2),dim(3));
nanMat(~nanIdx,:) = residual;
residual = reshape(nanMat,dim(1),dim(2),dim(3));

nanMat = NaN(dim(1)*dim(2),size(reg,2));
nanMat(~nanIdx,:) = rMat;
rMat = reshape(nanMat,dim(1),dim(2),[]);

% residual = sig - sum(permute(rMat,[1 2 4 3]).*reg,4);