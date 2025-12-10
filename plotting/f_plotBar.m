function [dataMean,dataSEM] = f_plotBar(data,varargin)

p = inputParser;
addParameter(p,'colors',[0, 0, 0]);
addParameter(p,'ylabel','');
addParameter(p,'legend',[]);
addParameter(p,'title','');
addParameter(p,'ylim',[]);
addParameter(p,'barVisible',1);
addParameter(p,'plotPoints',1);
addParameter(p,'norm',1);

parse(p,varargin{:});

plotLegend = p.Results.legend;
colors = p.Results.colors;

if iscell(plotLegend)
    plotLegend = [plotLegend; repmat({''},1,numel(plotLegend))];
    plotLegend = plotLegend(:);
end

N = numel(data);
dataMean = zeros(N,1);
dataSEM = zeros(N,1);
for i = 1:N
    dataMean(i) = mean(data{i});
    if p.Results.norm
        dataSEM(i) = std(data{i},0)/sqrt(numel(data{i}));
    else
        dataSEM(i) = std(data{i},0);
    end
end

hold on;
numZero = zeros(N,1);
if ~isempty(p.Results.ylim)
    ylim(p.Results.ylim);
    for i = 1:N
        numZero(i) = sum(data{i} < p.Results.ylim(1));
    end
end

cIdx = round(linspace(1,size(colors,1),numel(data)));

for i = 1:N
    if p.Results.barVisible
        b(i) = bar(i,dataMean(i));
    end
    if p.Results.plotPoints
        scatter(i*ones(numel(data{i}),1),data{i},70,'filled',XJitter='randn',XJitterWidth=0.3,MarkerFaceColor=colors(cIdx(i),:),MarkerFaceAlpha=0.5);
        if numZero(i)
            scatter(i*ones(numZero(i),1),zeros(numZero(i),1),100,'o',XJitter='randn',XJitterWidth=0.3,MarkerEdgeColor=colors(cIdx(i),:),MarkerFaceAlpha=0.5);
        end
    end

    if p.Results.barVisible
        b(i).FaceColor = 'flat';
        b(i).CData = [1 1 1];
        b(i).ShowBaseLine = 'off';
        b(i).EdgeColor = colors(cIdx(i),:);
        b(i).LineWidth = 3;
    end
end

er = errorbar(1:N,dataMean,dataSEM);

er.Color = [0 0 0];
er.LineStyle = 'none';
er.LineWidth = 3;

ax = gca;
ax.XAxis.Visible = 'off';

ylabel(p.Results.ylabel);
if ~isempty(p.Results.legend)
    legend(plotLegend);
end
title(p.Results.title);
set(gca,'FontSize',14);

end