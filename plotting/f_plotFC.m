function f_plotFC(varargin)

map = varargin{1};
diagVal = varargin{2};
map(diag(true(12,1))) = diagVal;

p = inputParser;
addParameter(p,'cmp',jet);
addParameter(p,'bounds',[0,1]);
addParameter(p,'clabel','');
addParameter(p,'title','');
addParameter(p,'cbar',1);

parse(p,varargin{3:end});

imagesc(map);
axis image off;
if p.Results.cbar
    c = colorbar;
    c.Label.String = p.Results.clabel;
end
ax = gca;
colormap(ax,p.Results.cmp);
clim(p.Results.bounds);
title(p.Results.title);
set(gca,'FontSize',14);

end