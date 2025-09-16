function [allen,overall,f] = f_hemLag_dT_allen(sig1,sig2,fs,lagwin,masks,toPlot)

%%
% sig1 = data.HbT;
% sig2 = data.rfp_HD;
% fs = 10;
% lagwin = [-5 5];
%

mask_all = sum(masks,3,'omitnan');
allen_labels = {'MOp','MOs','SSp-bfd','SSp-tr','SSp-ll','SSp-ul','SSp-un','VISpm','VISrl','VISam','VISa','VISp'};

dim = size(sig1);

sig1 = reshape(sig1,dim(1)*dim(2),dim(3));
sig2 = reshape(sig2,dim(1)*dim(2),dim(3));
nanIdx = mask_all==0;
sig1 = sig1(~nanIdx,:)';
sig2 = sig2(~nanIdx,:)';
sig1 = detrend(sig1,1);
sig2 = detrend(sig2,1);

maxlag = fs*max(abs(lagwin));
lagRange = -maxlag:maxlag;lagRange = lagRange >= lagwin(1)*fs & lagRange <= lagwin(2)*fs;

xcorrMat = f_xcorr(sig1,sig2,maxlag);

% calculate allen averages
n = size(masks,3);

maskIdx = sum(masks.*permute(1:n,[1 3 2]),3,'omitnan');
maskIdx = maskIdx(~(maskIdx==0));

allen = struct;
allen.xcorr = zeros(n,numel(lagRange));
allen.lag = (-maxlag:maxlag)/fs;

for idx = 1:n
    allen.xcorr(idx,:) = mean(xcorrMat(:,maskIdx==idx),2);
end

overall = struct;
overall.xcorr = mean(xcorrMat,2);
overall.lag = (-maxlag:maxlag)/fs;

%%
if toPlot
    f = figure('Position',[100 100 1000 600]);
    tl = f_tiledLayout(3,1);

    ax(1) = nexttile([2 1]);
    imagesc(allen.xcorr);
    c = colorbar;
    c.Label.String = 'r';
%     ax.YAxis.Label.Visible = 'on';
    ylabel(allen_labels');
    set(get(gca,'ylabel'),'Rotation',0,'VerticalAlignment','Middle','HorizontalAlignment','Right');
    set(ax(1),'YTick',[],'XTick',[]);
    colormap jet;
    clim([-1 1]);
    set(ax(1),'FontSize',16);
    
    ax(2) = nexttile;
    plot(overall.lag,overall.xcorr,'LineWidth',2);
    xlabel('Lag (s)');
    ylabel('r');
    set(ax(2),'FontSize',14);
    xlim(lagwin);
    box off;

else
    f = [];
end
