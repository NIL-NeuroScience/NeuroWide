function f_plotFC(varargin)

map = varargin{1};
diagVal = varargin{2};
map(diag(true(12,1))) = diagVal;

p = inputParser;
addParameter(p,'cmp',jet);
addParameter(p,'bounds',[0,1]);
addParameter(p,'clabel','');
addParameter(p,'title','');

parse(p,varargin{3:end});

imagesc(map);
axis image off;
c = colorbar;
colormap(p.Results.cmp);
clim(p.Results.bounds);
title(p.Results.title);
c.Label.String = p.Results.clabel;
set(gca,'FontSize',14);

end