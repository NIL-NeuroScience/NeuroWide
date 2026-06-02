function map(img, varargin)

p = inputParser;
addParameter(p, 'cmp', jet);
addParameter(p, 'bounds', prctile(img(:), [1,99]));
addParameter(p, 'clabel', '');
addParameter(p, 'title', '');

parse(p, varargin{:});

imAlpha = ~isnan(img);

hold on;

imagesc(img, AlphaData=imAlpha);
axis image off;
c = colorbar;
colormap(p.Results.cmp);
clim(p.Results.bounds);
title(p.Results.title);
c.Label.String = p.Results.clabel;
set(gca, FontSize=14);

end