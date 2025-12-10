function f_multiScatter(varargin)

X = varargin{1};
Y = varargin{2};

p = inputParser;
addParameter(p,'cmp',[0, 0, 0]);
addParameter(p,'ylabel','');
addParameter(p,'xlabel','');
addParameter(p,'title','');
addParameter(p,'ylim',[]);
addParameter(p,'xlim',[]);
addParameter(p,'alpha',1);
addParameter(p,'lineWidth',1);
addParameter(p,'plotAverage',1);
addParameter(p,'markerSize',50);
addParameter(p,'intercept',1);

parse(p,varargin{3:end});

N = numel(X);
cIdx = round(linspace(1,size(p.Results.cmp,1),N));

hold on;
for i = 1:N
    scatter(X{i},Y{i},p.Results.markerSize,'filled',MarkerFaceAlpha=p.Results.alpha,MarkerFaceColor=p.Results.cmp(cIdx(i),:));
end

if ~isempty(p.Results.xlim)
    xlim(p.Results.xlim);
    current_x = xlim;
else
    current_x = xlim;
    xlim(current_x);
end
if ~isempty(p.Results.ylim)
    ylim(p.Results.ylim);
    current_y = ylim;
else
    current_y = ylim;
    ylim(current_y);
end

% plot linear regression

if p.Results.intercept
    ends = zeros(N,2);
    for i = 1:N
        lm = fitlm(X{i},Y{i});
        lm = table2array(lm.Coefficients);
        ends(i,:) = current_x*lm(2,1)+lm(1,1);
        plot(current_x,ends(i,:),color=p.Results.cmp(cIdx(i),:),lineWidth=p.Results.lineWidth);
    end
    if p.Results.plotAverage
        plot(current_x,mean(ends),color=[0 0 0],lineWidth=2*p.Results.lineWidth);
    end
else
    slopes = zeros(N,1);
    for i = 1:N
        slopes(i) = X{i}(:) \ Y{i}(:);
        plot(current_x,current_x*slopes(i),color=p.Results.cmp(cIdx(i),:),lineWidth=p.Results.lineWidth);
    end
    if p.Results.plotAverage
        plot(current_x,current_x*mean(slopes),color=[0 0 0],lineWidth=2*p.Results.lineWidth);
    end
end

title(p.Results.title);
xlabel(p.Results.xlabel);
ylabel(p.Results.ylabel);
set(gca,'FontSize',14);

end