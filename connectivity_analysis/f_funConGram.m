function [fcGram,fIdx] = f_funConGram(sig,win)

dim = size(sig);

idx = 1:win(2):dim(1);
idx(idx-1+win(1) > dim(1)) = [];

fIdx = win(1)/2:win(2):dim(1)-win(1)/2;

fcGram = zeros(size(sig,2),size(sig,2),numel(idx));

for i = 1:numel(idx)
    fcGram(:,:,i) = corrcoef(sig(idx(i):idx(i)-1+win(1),:));
end