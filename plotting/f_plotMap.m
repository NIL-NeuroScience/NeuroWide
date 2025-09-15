function f_plotMap(varargin)

map = varargin{1};

p = inputParser;
addParameter(p,'cmp',jet);
addParameter(p,'bounds',prctile(map(:),[1,99]));
addParameter(p,'clabel','');
addParameter(p,'title','');

parse(p,varargin{2:end});

imAlpha = ~isnan(map);

hold on;

imagesc(map,AlphaData=imAlpha);
axis image off;
c = colorbar;
colormap(p.Results.cmp);
clim(p.Results.bounds);
title(p.Results.title);
c.Label.String = p.Results.clabel;
set(gca,'FontSize',14);

end