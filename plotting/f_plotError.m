function f_plotError(data,varargin)

p = inputParser;
addParameter(p,'color',[0, 0, 0]);
addParameter(p,'ylabel','');
addParameter(p,'title','');
addParameter(p,'ylim',[]);
addParameter(p,'lineWidth',2);
addParameter(p,'jitter',1);

parse(p,varargin{:});

[n,N] = size(data);

hold on;

dataMean = mean(data,1);
dataSEM = std(data,0,1)/sqrt(n);

er = errorbar(1:N,dataMean,dataSEM);
er.Color = p.Results.color;
er.LineWidth = p.Results.lineWidth;

for i = 1:N
    if p.Results.jitter
        s(1) = scatter(i*ones(n,1),data(:,i),100,'filled',MarkerFaceColor=p.Results.color,XJitter='randn',XJitterWidth=0.3,MarkerFaceAlpha=0.5);
    else
        s(1) = scatter(i*ones(n,1),data(:,i),100,'filled',MarkerFaceColor=p.Results.color,MarkerFaceAlpha=0.5);
    end
end
ax = gca;
ax.XAxis.Visible = 'off';
ylabel(p.Results.ylabel);

title(p.Results.title);
set(gca,'FontSize',14);

end