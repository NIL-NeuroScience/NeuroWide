function [allen,overall,f] = f_hemSpectra_allen(sig,fs,fpass,tapers,masks,toPlot,fLim)

%%
% sig = data.rfp_HD;
% fs = 10;
% fpass = [0 5];
% tapers = [5 9];
%

mask_all = sum(masks,3,'omitnan');
allen_labels = {'MOp','MOs','SSp-bfd','SSp-tr','SSp-ll','SSp-ul','SSp-un','VISpm','VISrl','VISam','VISa','VISp'};

dim = size(sig);
sig = reshape(sig,[dim(1)*dim(2) dim(3)]);
nanIdx = mask_all==0;
sig = sig(~nanIdx,:);
sig = sig';

params.Fs = fs;
params.fpass = fpass;
params.trialave = 0;
params.tapers = tapers;

[spectraMat,fr] = mtspectrumc(sig,params);

% calculate allen averages
n = size(masks,3);

maskIdx = sum(masks.*permute(1:n,[1 3 2]),3,'omitnan');
maskIdx = maskIdx(~(maskIdx==0));

allen = struct;
allen.spectra = zeros(n,numel(fr));
allen.fr = fr;

for idx = 1:n
    allen.spectra(idx,:) = mean(spectraMat(:,maskIdx==idx),2);
end

overall = struct;
overall.spectra = mean(spectraMat,2);
overall.fr = fr;

%%
if toPlot
    f = figure('Position',[100 100 1000 600]);
    tl = f_tiledLayout(3,1);

    ax(1) = nexttile([2 1]);
    imagesc(log10(allen.spectra));
    c = colorbar;
    c.Label.String = 'Log_1_0(Power)';
%     ax.YAxis.Label.Visible = 'on';
    ylabel(allen_labels');
    set(get(gca,'ylabel'),'Rotation',0,'VerticalAlignment','Middle','HorizontalAlignment','Right');
    set(ax(1),'YTick',[],'XTick',0:numel(fr)/10:numel(fr),'XTickLabel',0:fs/20:fs/2);
    colormap cmpinf;
    clim(prctile(log10(allen.spectra(:)),[1 99]));
    set(ax(1),'FontSize',16);

    if nargin > 6
        [~,tmpfLim] = min(fr'-fLim,[],1,'ComparisonMethod','abs');
        xlim(tmpfLim);
        tmp = allen.spectra(:,tmpfLim(1):tmpfLim(2));
        clim(prctile(log10(tmp(:)),[1 99]));
    end
    
    ax(2) = nexttile;
    plot(fr,log10(overall.spectra),'LineWidth',2);
    xlabel('Frequency (Hz)');
    ylabel('Log_1_0(Power)');
    set(ax(2),'FontSize',14);
    box off; 
    if nargin > 6
        xlim(fLim);
    end

else
    f = [];
end


