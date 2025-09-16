function mask = f_combineParcells(parcells)

mask = nan(size(parcells(:,:,:,1)));

N = size(parcells,3);

for idx = 1:N
    tmp = nan(size(mask(:,:,1)));
    tmp(parcells(:,:,idx,1) == 1) = 1;
    tmp(parcells(:,:,idx,2) == 1) = 1;
    mask(:,:,idx) = tmp;
end
