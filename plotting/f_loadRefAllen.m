function f_loadRefAllen()
%%
path = fullfile(f_path,'plotting/refAllen.mat');

tmp = load(path);

%%
refBM = tmp.refBM;
refParcellation = tmp.refParcellation;
refMasks = refParcellation.Masks;


