function f_plotBox(data,varargin)

p = inputParser;
addParameter(p,'cmp',[0, 0, 0]);
addParameter(p,'ylabel','');
addParameter(p,'title','');
addParameter(p,'ylim',[]);
addParameter(p,'lineWidth',1);

parse(p,varargin{:});

[n,N] = size(data);

hold on;

cIdx = round(linspace(1,size(p.Results.cmp,1),N));

for i = 1:N
    s(1) = scatter(i*ones(n,1),data(:,i),50,'filled',MarkerFaceColor=p.Results.cmp(cIdx(i),:),XJitter='randn',XJitterWidth=0.3,MarkerFaceAlpha=0.5);
    b(i) = boxchart(i*ones(n,1),data(:,i),MarkerStyle='none');
    b(i).BoxFaceColor = [1,1,1];
    b(i).BoxEdgeColor = p.Results.cmp(cIdx(i),:);
    b(i).LineWidth = p.Results.lineWidth;
end

ax = gca;
ax.XAxis.Visible = 'off';
ylabel(p.Results.ylabel);

title(p.Results.title);
set(gca,'FontSize',14);

end