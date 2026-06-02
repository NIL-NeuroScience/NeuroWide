function [img,fig_img] = allenMap(data, varargin)

p = inputParser;
addParameter(p, 'mask', []);
addParameter(p, 'parcellation', []);
addParameter(p, 'cmp', []);
addParameter(p, 'title', []);
addParameter(p, 'cLabel', []);
addParameter(p, 'clim', []);
addParameter(p, 'cbar', 1);
addParameter(p, 'plot', 1);
addParameter(p, 'ylabel', []);

parse(p, varargin{:});

parcellation = p.Results.parcellation;
mask = p.Results.mask;
cmp = p.Results.cmp;
Title = p.Results.title;
cLabel = p.Results.cLabel;
cRange = p.Results.clim;

if isempty(parcellation) && isempty(mask)
    [mask,parcellation] = NeuroWide.plot.loadRefAllen();
elseif isempty(parcellation)
    [~,parcellation] = NeuroWide.plot.loadRefAllen();
end

masks = sum(parcellation.Masks, 4);

if ~isempty(mask)
    mask(isnan(mask)) = 0;
    masks = masks .* mask;
end

masks = logical(masks);
img = NaN(size(masks, [1,2]));

for i = 1:size(masks, 3)
    img(masks(:,:,i)) = data(i);
end

imAlpha = ~isnan(img);

y_crop = find(sum(imAlpha, 2));
y_crop = y_crop(1):y_crop(end);
x_crop = find(sum(imAlpha, 1));
x_crop = x_crop(1):x_crop(end);

imAlpha = imAlpha(y_crop, x_crop);

if p.Results.plot
    fig_img = imagesc(img(y_crop, x_crop), AlphaData=imAlpha);
    axis image off;
    if ~isempty(cmp)
        ax = gca;
        colormap(ax, cmp);
    end
    
    if p.Results.cbar
        c = colorbar;
        c.Label.String = cLabel;
    end
    title(Title);
    set(gca, FontSize=14);
    if ~isempty(cRange)
        clim(cRange);
    end

    if ~isempty(p.Results.ylabel)
        ylabel(p.Results.ylabel);
        ax = gca;
        ax.YLabel.Visible = 'on';
    end
end