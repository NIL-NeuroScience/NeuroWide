function f_plotAllenRegion(region,varargin)

p = inputParser;
addParameter(p,'mask',[]);
addParameter(p,'parcellation',[]);
addParameter(p,'color',[0, 0, 0]);
addParameter(p,'linewidth',1);

parse(p,varargin{2:end});

if numel(varargin) > 0
    side = varargin{1};
else
    side = [];
end

parcellation = p.Results.parcellation;
mask = p.Results.mask;
color = p.Results.color;
linewidth = p.Results.linewidth;

allen_path = fullfile(f_path,'plotting/refAllen.mat');

if isempty(parcellation) && isempty(mask)
    parcellation = load(allen_path);
    mask = parcellation.refBM;
    parcellation = parcellation.refParcellation;
elseif isempty(parcellation)
    parcellation = load(allen_path);
    parcellation = parcellation.refParcellation;
end

if isempty(side)
    masks = sum(parcellation.Masks(:,:,region,:),4);
else
    masks = parcellation.Masks(:,:,region,side);
end

if ~isempty(mask)
    mask(isnan(mask)) = 0;
    masks = masks.*mask;
end

bound = bwboundaries(masks);
N = numel(bound);

hold on;
for i = 1:N
    plot(bound{i}(:,2),bound{i}(:,1),Color=color,LineWidth=linewidth);
end

set(gca,'YDir','reverse');