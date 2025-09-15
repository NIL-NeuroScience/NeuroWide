function f_plotAllenMap(varargin)

p = inputParser;
addParameter(p,'mask',[]);
addParameter(p,'parcellation',[]);
addParameter(p,'cmp',[]);
addParameter(p,'title',[]);
addParameter(p,'cLabel',[]);
addParameter(p,'cRange',[]);

parse(p,varargin{2:end});

data = varargin{1};

parcellation = p.Results.parcellation;
mask = p.Results.mask;
cmp = p.Results.cmp;
Title = p.Results.title;
cLabel = p.Results.cLabel;
cRange = p.Results.cRange;

allen_path = fullfile(f_path,'plotting/refAllen.mat');

if isempty(parcellation) && isempty(mask)
    parcellation = load(allen_path);
    mask = parcellation.refBM;
    parcellation = parcellation.refParcellation;
elseif isempty(parcellation)
    parcellation = load(allen_path);
    parcellation = parcellation.refParcellation;
end

masks = sum(parcellation.Masks,4);

if ~isempty(mask)
    mask(isnan(mask)) = 0;
    masks = masks.*mask;
end

masks = logical(masks);
img = NaN(size(masks,[1,2]));

for i = 1:size(masks,3)
    img(masks(:,:,i)) = data(i);
end

imAlpha = ~isnan(img);

imagesc(img,AlphaData=imAlpha);
axis image off;
if ~isempty(cmp)
    colormap(cmp);
end
c = colorbar;
c.Label.String = cLabel;
title(Title);
set(gca,'FontSize',14);
if ~isempty(cRange)
    clim(cRange);
end