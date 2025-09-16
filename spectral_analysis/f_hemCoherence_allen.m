function [allen,overall,f] = f_hemCoherence_allen(sig1,sig2,fs,masks,tapers,toPlot,fLim)

%%
% sig1 = data.HbT;
% sig2 = data.gfp_HD;
% fs = 10;
% tapers = [5 9];
%

mask_all = sum(masks,3,'omitnan');
allen_labels = {'MOp','MOs','SSp-bfd','SSp-tr','SSp-ll','SSp-ul','SSp-un','VISpm','VISrl','VISam','VISa','VISp'};

dim = size(sig1);
nanIdx = mask_all==0;
sig1 = reshape(sig1,dim(1)*dim(2),[]);
sig2 = reshape(sig2,dim(1)*dim(2),[]);
sig1 = sig1(~nanIdx,:)';
sig2 = sig2(~nanIdx,:)';

params = struct;
params.Fs = fs;
params.tapers = tapers;
params.trialave = 0;

[C,phi,~,~,~,fr] = coherencyc(sig1,sig2,params);

n = size(masks,3);

maskIdx = sum(masks.*permute(1:n,[1 3 2]),3,'omitnan');
maskIdx = maskIdx(~(maskIdx==0));

allen = struct;
allen.C = zeros(n,numel(fr));
allen.phi = allen.C;
allen.fr = fr;

for idx = 1:n
    allen.C(idx,:) = mean(C(:,maskIdx==idx),2);
    allen.phi(idx,:) = mean(phi(:,maskIdx==idx),2);
end

overall = struct;
overall.C = mean(C,2);
overall.phi = mean(phi,2);
overall.fr = fr;

%%
if toPlot
    f = figure('Position',[100 100 1000 800]);
    tl = f_tiledLayout(2,1);

    ax(1) = nexttile;
    imagesc(allen.C);
    c = colorbar('XTick',[0 1]);
    c.Label.String = 'coherence';
%     ax.YAxis.Label.Visible = 'on';
    ylabel(allen_labels');
%     xlabel('Frequency (Hz)');
    set(get(ax(1),'ylabel'),'Rotation',0,'VerticalAlignment','Middle','HorizontalAlignment','Right');
    set(ax(1),'YTick',[],'XTick',0:numel(fr)/10:numel(fr),'XTickLabel',0:fs/20:fs/2);
    colormap cmpinf;
    clim([0 1]);
    set(gca,'FontSize',16);

    if nargin > 6
        [~,fLim] = min(fr'-fLim,[],1,'ComparisonMethod','abs');
        xlim(fLim);
    end

    ax(2) = nexttile;
    imagesc(allen.phi);
    c = colorbar('XTick',[-3.14 3.14]);
    c.Label.String = 'rad';
%     ax.YAxis.Label.Visible = 'on';
    ylabel(allen_labels');
    xlabel('Frequency (Hz)');
    set(get(ax(2),'ylabel'),'Rotation',0,'VerticalAlignment','Middle','HorizontalAlignment','Right');
    set(ax(2),'YTick',[],'XTick',0:numel(fr)/10:numel(fr),'XTickLabel',0:fs/20:fs/2);
    colormap jet;
    clim([-pi pi]);
    set(gca,'FontSize',16);

    if nargin > 6
        xlim(fLim);
    end

    colormap(ax(1),cmpinf);
else
    f = [];
end





