function f_plotLineError(varargin)

x = varargin{1}(:);
y = varargin{2}(:);
error = varargin{3}(:);

p = inputParser;
addParameter(p,'color',[]);
addParameter(p,'lineWidth',2);
addParameter(p,'log',0);

parse(p,varargin{4:end});

if isempty(p.Results.color)
    color = get(groot,'defaultAxesColorOrder');
    color = color(1,:);
else
    color = p.Results.color;
end

if p.Results.log
    zeroIdx = x == 0;
    x(zeroIdx) = [];
    y(zeroIdx) = [];
    error(zeroIdx) = [];
    set(gca,'YScale','log','XScale','log');
end

hold on;
fill([x;flipud(x)],[(y+error);flipud(y-error)],color,FaceAlpha=0.3,EdgeColor='none');
plot(x,y,Color=color,LineWidth=p.Results.lineWidth);

end