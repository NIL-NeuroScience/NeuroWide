function parcellation = f_parcellate(sig,masks)

masks = reshape(masks,[],size(masks,3));

masks(isnan(masks)) = 0;

dim = size(sig);

sig = reshape(sig,[],dim(3))';
sig(isnan(sig)) = 0;

parcellation = sig*masks;
parcellation = parcellation./sum(masks,1);

end