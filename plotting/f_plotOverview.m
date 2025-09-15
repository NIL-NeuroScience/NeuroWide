function f = f_plotOverview(Signals,GRAB)

pS = struct;
pS.fontSize = 14;
pS.xTick = 0:100:Signals.t(end);
pS.lW = 1;
pS.rfp_exist = isfield(Signals,'rfp_HD');
pS.gfp_exist = isfield(Signals,'gfp_HD');
pS.xlim = [0,Signals.t(end)];

pS.title = ['\color[rgb]{0.8196, 0.0784, 0.0941}Ca^2^+   \color[rgb]{0.4471, 0.7529, 0.0314}' GRAB '   \color[rgb]{0.2039 0.2157 0.5922}HbT   \color[rgb]{0.7 0 0.7}Pupil Diameter   \color[rgb]{0 0.7 0.7}Whisking   \color[rgb]{0 0 0}Movement   \color[rgb]{0.4 0.4 0.4}EMG'];

% find max and min

shifts = [0,10,20];
regions = [2,5,3,12];

all_plot = [];

if pS.rfp_exist
    rfp_plot = Signals.rfp_HD(:,regions);
    all_plot = rfp_plot-shifts(1);
end
if pS.gfp_exist
    gfp_plot = Signals.gfp_HD(:,regions);
    all_plot = [all_plot gfp_plot-shifts(2)];
end
HbT_plot = Signals.HbT(:,regions);
all_plot = [all_plot HbT_plot-shifts(3)];

pS.limit = prctile(all_plot,[0.1,99.9]);
pS.limit = [-1;1].*max(pS.limit.*[-1;1],[],2);

% plot

regions = {'MOs','SSp-ll','SSp-bfd','VISp'};

f = figure('Position',[0 0 1200 900]);

tl = tiledlayout(11,1);
tl.TileSpacing = 'compact';
tl.Padding = 'compact';
title(tl,pS.title,'interpreter','tex','FontSize',pS.fontSize,'FontWeight','bold');

for i = 1:4
    ax(i) = nexttile([2,1]);hold on;
    plot(Signals.t,all_plot(:,i),Color=c_Ca,LineWidth=pS.lW);
    plot(Signals.t,all_plot(:,4+i),Color=c_GRAB,LineWidth=pS.lW);
    plot(Signals.t,all_plot(:,8+i),Color=c_HbT,LineWidth=pS.lW);
    if i == 1
        plot([Signals.t(end) Signals.t(end)],[pS.limit(2) pS.limit(2)-10],'-k',LineWidth=3);
    end
    box off;
    ax(i).XAxis.Visible = 'off';
    ax(i).YAxis.Visible = 'off';
    ax(i).YAxis.Label.Visible = 'on';
    ylabel(regions{i});
    ylim(pS.limit);
    xlim(pS.xlim);
end

ax(5) = nexttile([3,1]);hold on;

if isfield(Signals,'pupil') && ~isempty(Signals.pupil)
    plot(Signals.t,Signals.pupil*0.4+0.6,color=c_pupil,LineWidth=pS.lW);
    plot([0,0],[0.6,1],'-k',lineWidth=3);
end
if isfield(Signals,'whisker_pad') && ~isempty(Signals.whisker_pad)
    plot(Signals.t,rescale(Signals.whisker_pad,0.3,0.6),color=[0,0.7,0.7],LineWidth=pS.lW);
end
if isfield(Signals,'Acc') && ~isempty(Signals.Acc)
    plot(Signals.t,rescale(Signals.Acc,0,0.3),color=[0,0,0],LineWidth=pS.lW);
end
ax(5).YAxis.Visible = 'off';
xlim(pS.xlim);
xlabel('Time (s)');
box off;
ax(5).YAxis.Label.Visible = 'on';
ylabel('Behavior');

set(ax(1:5),FontSize=pS.fontSize);

end