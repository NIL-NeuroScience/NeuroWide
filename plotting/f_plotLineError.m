function f_plotLineError(x,y,error,varargin)

x = x(:);
y = y(:);
[h,w] = size(error);

if h == 1 & w > 1
    error = error(:);
end

p = inputParser;
addParameter(p,'color',[]);
addParameter(p,'ylabel',[]);
addParameter(p,'xlim',[]);
addParameter(p,'lineWidth',2);
addParameter(p,'log',0);

parse(p,varargin{:});

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
if numel(error) == numel(x)
    fill([x;flipud(x)],[(y+error);flipud(y-error)],color,FaceAlpha=0.3,EdgeColor='none');
else
    fill([x;flipud(x)],[error(:,2);flipud(error(:,1))],color,FaceAlpha=0.3,EdgeColor='none');
end
plot(x,y,Color=color,LineWidth=p.Results.lineWidth);

if ~isempty(p.Results.ylabel)
    ylabel(p.Results.ylabel);
end

if ~isempty(p.Results.xlim)
    xlim(p.Results.xlim);
end

end