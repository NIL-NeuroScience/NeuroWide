%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%           Single photon data postprocessing and analysis
%
% Secondary processing script for widefield fluorescence and intrinsic
% imaging data. Creates brain and vessel masks, estimates pupil diameter,
% whisking, and movement, and performs preliminary analysis between the
% optical signals.
% 
% INPUTS: f_ROIprocessing(Mouse,Date,Runs,_)
%   Mouse: mouse ID
%   Date: session date
%   Runs: list of all runs to be analyzed. Leave empty to analyze all 
%       processed runs
% 
% EXTRA PARAMETERS:
%   behCam: true if behavior was recorded (default = true)
%   newROIs: true to define new ROIs (default = false)
%   extraProcessing: true to perform extra analysis (default = true)
%   GRAB: description of gfp channel (default = 'none')
%   load_dir: directory to load from (default = 'bcraus/HRF/1P')
%   save_dir: directory to save to (default = 'bcraus/HRF/1P')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f_ROIprocessing(Mouse,Date,Runs,varargin)
%% parse inputs
p = inputParser;
addParameter(p,'behCam',true);
addParameter(p,'newROIs',false);
addParameter(p,'extraProcessing',true);
addParameter(p,'GRAB','none');
addParameter(p,'load_dir','bcraus/HRF/1P');
addParameter(p,'save_dir','bcraus/HRF/1P');
addParameter(p,'plotVisible',false);

parse(p,varargin{:});

%% organize files

files = struct;
files.load_dir = fullfile('/projectnb/devorlab',p.Results.load_dir);
files.save_dir = fullfile('/projectnb/devorlab',p.Results.save_dir);
files.load_dir = fullfile(files.load_dir,Date,Mouse);
files.save_dir = fullfile(files.save_dir,Date,Mouse,'DataAnalysis');
files.dataIn = fullfile(files.load_dir,'dataIn.mat');
files.triggers = fullfile(files.load_dir,'Triggers');
files.images = fullfile(files.save_dir,'Images');
files.processed = fullfile(files.load_dir,'processed');
files.camera = fullfile(files.load_dir,'camera');

[~,~,~,] = mkdir(files.save_dir);
[~,~,~] = mkdir(files.images);

dataIn = load(files.dataIn);
dataIn = dataIn.dataIn;

%% find runs
list_runs = dir(files.processed);
list_runs = {list_runs(3:end).name};
list_runs = cellfun(@(f) str2double(f(4:7)),list_runs);

if ~isempty(Runs)
    list_runs = intersect(list_runs,Runs);
end

%% create brain and behavior masks

ROI_exist = exist(fullfile(files.save_dir,'ROIs.mat'),'file');

dataIn_idx = [dataIn.runnum]==list_runs(1);

if p.Results.newROIs || ~ROI_exist
    [brain_mask,hem] = f_createROIs(dataIn(dataIn_idx).template(:,:,1));
    
    % create vessel mask
    vessel_good = false;
    sensitivity = 15;

    LEDs = [dataIn(dataIn_idx).led(:).type];
    idx_525 = find(LEDs==525);
    f1 = figure(Position=[400,300,700,600]);axis image off;hold on;colormap gray;set(gca,'YDir','reverse');
    imagesc(dataIn(dataIn_idx).template(:,:,idx_525));
    clim([0,65535]);

    f = figure(Position=[400,300,700,600]);axis image off;hold on;colormap gray;set(gca,'YDir','reverse');
    while ~vessel_good
        vessel_mask = f_maskVessel(dataIn(dataIn_idx).template(:,:,idx_525),sensitivity,brain_mask);
        imagesc(dataIn(dataIn_idx).template(:,:,idx_525));
        clim([0,65535]);
        imagesc(vessel_mask,AlphaData=0.5*isnan(vessel_mask));
        
        vessel_good = input('Good?');
        if ~vessel_good
            sensitivity = input('New sensitivity? (default = 15)');
        end
    end
    close(f1);
    close(f);
    
    save(fullfile(files.save_dir,'ROIs.mat'),'brain_mask','vessel_mask','hem');

    if p.Results.behCam
        cam_folder = fullfile(files.camera,sprintf('Run%02i',list_runs(1)));
        cam_images = dir(cam_folder);cam_images(1:2) = [];
        cam_images = fullfile({cam_images.folder},{cam_images.name});

        imgPath = cam_images{round((numel(cam_images))/2)};
        behaviorROIs = f_drawBehaviorROIs(imgPath);
        clear cam_folder cam_images imgPath;
        
        save(fullfile(files.save_dir,'ROIs.mat'),'behaviorROIs','-append');        
    end
else
    load(fullfile(files.save_dir,'ROIs.mat'));
end

%% create allen masks

allen_exist = who('-file',fullfile(files.save_dir,'ROIs.mat'));
allen_exist = ismember('parcellation',allen_exist);

if p.Results.newROIs || ~allen_exist
    LEDs = [dataIn(dataIn_idx).led(:).type];
    rfp_idx = find(LEDs==565);
    
    if isempty(rfp_idx)
        rfp_idx = 1;
    end
    
    importSettings.mp = 1;
    importSettings.cores= 4;
    importSettings.doNorm = 0;
    importSettings.doCat=0;
    importSettings.nImport = 0;
    importSettings.nChannels = numel(LEDs);
    
    rawData = f_AndorDATImporter(fullfile(files.load_dir,'onephoton',sprintf('Run%02i',list_runs(1))),importSettings);
    rawData = std(rawData(:,:,:,rfp_idx),0,3);

    rotation = dataIn(1).rotation;
    
    parcellation = f_AllenAtlas(imrotate(rawData,rotation),brain_mask,hem);
    
    save(fullfile(files.save_dir,'ROIs.mat'),'parcellation','rotation','-append');
end

%%
for Run = list_runs
    %%
    disp(['Analyzing ' Date ', ' Mouse ', Run ' num2str(Run)]);

    metadata = struct;

    % load triggers and data
    [~,~,~] = mkdir(fullfile(files.save_dir,sprintf('Run%02i',Run)));
    settings = load(fullfile(files.load_dir,'Triggers',sprintf('Run%03i.mat',Run)));
    digitalInput = settings.digitalInput;
    settings = settings.settings;

    % load h5 files
    files.h5 = fullfile(files.processed,sprintf('run%04i.h5',Run));

    LEDs = [dataIn(dataIn_idx).led(:).type];
    rfp_exist = sum(LEDs==565);
    gfp_exist = sum(LEDs==470);

    masks = sum(parcellation.Masks,4);

    Signals = struct;
    if gfp_exist
        gfp_HD = 100*h5read(files.h5,'/gfp/normHD');
        gfp = 100*h5read(files.h5,'/gfp/norm');
        Signals.gfp_HD = f_parcellate(gfp_HD,masks.*vessel_mask);
        Signals.gfp = f_parcellate(gfp,masks.*vessel_mask);
    end
    if rfp_exist
        rfp_HD = 100*h5read(files.h5,'/rfp/normHD');
        rfp = 100*h5read(files.h5,'/rfp/norm');
        Signals.rfp_HD = f_parcellate(rfp_HD,masks.*vessel_mask);
        Signals.rfp = f_parcellate(rfp,masks.*vessel_mask);
    end
    HbO = 1e6*h5read(files.h5,'/hemodynamics/HbO');
    HbR = 1e6*h5read(files.h5,'/hemodynamics/Hb');
    HbT = HbO+HbR;
    Signals.HbO = f_parcellate(HbO,masks);
    Signals.HbR = f_parcellate(HbR,masks);
    Signals.HbT = f_parcellate(HbT,masks);

    %% process behavior

    Signals.Acc = f_alignAcc(digitalInput.SonaTrigger,digitalInput.Accelerometer);
    sampling_win = floor(settings.DAQFrequency/settings.fs);
    Signals.Acc = movstd(Signals.Acc,sampling_win,0);
    Signals.Acc = Signals.Acc(round(sampling_win/2):sampling_win:end);

    Signals.t = 1/settings.fs:1/settings.fs:settings.totalTime;

    if p.Results.behCam
    % behavior extraction
        cam_folder = fullfile(files.camera,sprintf('Run%02i',Run));
        cam_images = dir(cam_folder);cam_images(1:2) = [];
        cam_images = fullfile({cam_images.folder},{cam_images.name});
    
        % pupil
        Signals.pupil = f_processPupil(cam_images,behaviorROIs);
    
        % whisker
        [Signals.whisker_long,Signals.whisker_pad] = f_processWhisking(cam_images,behaviorROIs);
    
        beh_padding = ones(1,(settings.nCycles-size(Signals.pupil,2)));
        Signals.pupil = [Signals.pupil Signals.pupil(end)*beh_padding];
        Signals.whisker_long = [Signals.whisker_long Signals.whisker_long(end)*beh_padding];
        Signals.whisker_pad = [Signals.whisker_pad Signals.whisker_pad(end)*beh_padding];
    end

    % save Signals
    save(fullfile(files.save_dir,sprintf('Run%02i',Run),'Signals.mat'),'Signals');
    
    %% create images
    
    iRun_length = 600*settings.fs;
    run_bins = floor(numel(Signals.t)/iRun_length);
    run_idx = {};
    for i = 1:run_bins
        if i == run_bins
            run_idx{i} = (i-1)*iRun_length+1:numel(Signals.t);
        else
            run_idx{i} = (i-1)*iRun_length+1:i*iRun_length;
        end
    end
    
    % Overview plots
    
    if ~p.Results.plotVisible
        set(0,'DefaultFigureVisible','off');
    end

    f = f_plotOverview(Signals,p.Results.GRAB);
    exportgraphics(f,fullfile(files.images,sprintf('plotOverview_Run%02i.png',Run)),'Resolution',300,'BackgroundColor',[1 1 1]);
    savefig(f,fullfile(files.images,sprintf('plotOverview_Run%02i.fig',Run))); close(f);
        
        % Correlation images
    for i = 1:run_bins
        iR_metadata = struct;
        if gfp_exist
            [iR_metadata.rMaps.gfpHbT,f] = f_corr(brain_mask.*gfp_HD(:,:,run_idx{i}),HbT(:,:,run_idx{i}),3,plot=true);clim([-1 1]);axis image off;
            title('GFP vs. HbT');
            img_name = sprintf('gfpHbTCorr_Run%02i_iRun%02i',Run,i);
            exportgraphics(f,fullfile(files.images,[img_name '.png']),'Resolution',300,'BackgroundColor',[1 1 1]);
            savefig(f,fullfile(files.images,[img_name '.fig'])); close(f);
        end
    
        if rfp_exist
            if gfp_exist
                [iR_metadata.rMaps.rfpgfp,f] = f_corr(brain_mask.*gfp_HD(:,:,run_idx{i}),rfp_HD(:,:,run_idx{i}),3,plot=true);clim([-1 1]);axis image off;
                title('GFP vs. RFP');
                img_name = sprintf('gfprfpCorr_Run%02i_iRun%02i',Run,i);
                exportgraphics(f,fullfile(files.images,[img_name '.png']),'Resolution',300,'BackgroundColor',[1 1 1]);
                savefig(f,fullfile(files.images,[img_name '.fig'])); close(f);
            end
            [iR_metadata.rMaps.rfpHbT,f] = f_corr(brain_mask.*rfp_HD(:,:,run_idx{i}),HbT(:,:,run_idx{i}),3,plot=true);clim([-1 1]);axis image off;
            title('RFP vs. HbT');
            img_name = sprintf('rfpHbTCorr_Run%02i_iRun%02i',Run,i);
            exportgraphics(f,fullfile(files.images,[img_name '.png']),'Resolution',300,'BackgroundColor',[1 1 1]);
            savefig(f,fullfile(files.images,[img_name '.fig'])); close(f);
        end
        
        %
        if p.Results.extraProcessing        
            % calculate coherence
            if gfp_exist
                [allen,overall,f] = f_hemCoherence_allen(gfp_HD(:,:,run_idx{i}),HbT(:,:,run_idx{i}),settings.fs,masks.*vessel_mask,[5 9],1,[0 0.5]);
                iR_metadata.coherence.fr = allen.fr;
                iR_metadata.coherence.allen.C.gfp_HbT = allen.C;
                iR_metadata.coherence.allen.phi.gfp_HbT = allen.phi;
                iR_metadata.coherence.overall.C.gfp_HbT = overall.C;
                iR_metadata.coherence.overall.phi.gfp_HbT = overall.phi;
                img_name = sprintf('coh_gfpHbT_Run%02i_iRun%02i',Run,i);
                exportgraphics(f,fullfile(files.images,[img_name '.png']),'Resolution',300,'BackgroundColor',[1 1 1]);
                savefig(f,fullfile(files.images,[img_name '.fig'])); close(f);
                if rfp_exist
                    [allen,overall,f] = f_hemCoherence_allen(gfp_HD(:,:,run_idx{i}),rfp_HD(:,:,run_idx{i}),settings.fs,masks.*vessel_mask,[5 9],1,[0 0.5]);
                    iR_metadata.coherence.fr = allen.fr;
                    iR_metadata.coherence.allen.C.rfp_gfp = allen.C;
                    iR_metadata.coherence.allen.phi.rfp_gfp = allen.phi;
                    iR_metadata.coherence.overall.C.rfp_gfp = overall.C;
                    iR_metadata.coherence.overall.phi.rfp_gfp = overall.phi;
                    img_name = sprintf('coh_rfpgfp_Run%02i_iRun%02i',Run,i);
                    exportgraphics(f,fullfile(files.images,[img_name '.png']),'Resolution',300,'BackgroundColor',[1 1 1]);
                    savefig(f,fullfile(files.images,[img_name '.fig'])); close(f);
                end
            end
            if rfp_exist
                [allen,overall,f] = f_hemCoherence_allen(rfp_HD(:,:,run_idx{i}),HbT(:,:,run_idx{i}),settings.fs,masks.*vessel_mask,[5 9],1,[0 0.5]);
                iR_metadata.coherence.fr = allen.fr;
                iR_metadata.coherence.allen.C.rfp_HbT = allen.C;
                iR_metadata.coherence.allen.phi.rfp_HbT = allen.phi;
                iR_metadata.coherence.overall.C.rfp_HbT = overall.C;
                iR_metadata.coherence.overall.phi.rfp_HbT = overall.phi;
                img_name = sprintf('coh_rfpHbT_Run%02i_iRun%02i',Run,i);
                exportgraphics(f,fullfile(files.images,[img_name '.png']),'Resolution',300,'BackgroundColor',[1 1 1]);
                savefig(f,fullfile(files.images,[img_name '.fig'])); close(f);
            end
    
            % calculate spectra
            if gfp_exist
                [allen,overall,f] = f_hemSpectra_allen(gfp_HD(:,:,run_idx{i}),settings.fs,[0 settings.fs/2],[5 9],masks.*vessel_mask,1,[0 0.5]);
                iR_metadata.spectra.fr = allen.fr;
                iR_metadata.spectra.allen.gfp = allen.spectra;
                iR_metadata.spectra.overall.gfp = overall.spectra;
                img_name = sprintf('spectra_gfp_Run%02i_iRun%02i',Run,i);
                exportgraphics(f,fullfile(files.images,[img_name '.png']),'Resolution',300,'BackgroundColor',[1 1 1]);
                savefig(f,fullfile(files.images,[img_name '.fig'])); close(f);
            end
            if rfp_exist
                [allen,overall,f] = f_hemSpectra_allen(rfp_HD(:,:,run_idx{i}),settings.fs,[0 settings.fs/2],[5 9],masks.*vessel_mask,1,[0 5]);
                iR_metadata.spectra.fr = allen.fr;
                iR_metadata.spectra.allen.rfp = allen.spectra;
                iR_metadata.spectra.overall.rfp = overall.spectra;
                img_name = sprintf('spectra_rfp_Run%02i_iRun%02i',Run,i);
                exportgraphics(f,fullfile(files.images,[img_name '.png']),'Resolution',300,'BackgroundColor',[1 1 1]);
                savefig(f,fullfile(files.images,[img_name '.fig'])); close(f);
            end
            [allen,overall,f] = f_hemSpectra_allen(HbT(:,:,run_idx{i}),settings.fs,[0 settings.fs/2],[5 9],masks,1,[0 0.5]);
            iR_metadata.spectra.fr = allen.fr;
            iR_metadata.spectra.allen.HbT = allen.spectra;
            iR_metadata.spectra.overall.HbT = overall.spectra;
            img_name = sprintf('spectra_HbT_Run%02i_iRun%02i',Run,i);
            exportgraphics(f,fullfile(files.images,[img_name '.png']),'Resolution',300,'BackgroundColor',[1 1 1]);
            savefig(f,fullfile(files.images,[img_name '.fig'])); close(f);
    
            % calculate xcorr
            if gfp_exist
                [allen,overall,f] = f_hemLag_dT_allen(HbT(:,:,run_idx{i}),gfp_HD(:,:,run_idx{i}),settings.fs,[-5 5],masks.*vessel_mask,1);
                iR_metadata.xcorr.lag = allen.lag;
                iR_metadata.xcorr.allen.HbT_gfp = allen.xcorr;
                iR_metadata.xcorr.overall.HbT_gfp = overall.xcorr;
                img_name = sprintf('xcorr_HbTgfp_Run%02i_iRun%02i',Run,i);
                exportgraphics(f,fullfile(files.images,[img_name '.png']),'Resolution',300,'BackgroundColor',[1 1 1]);
                savefig(f,fullfile(files.images,[img_name '.fig'])); close(f);
                if rfp_exist
                    [allen,overall,f] = f_hemLag_dT_allen(gfp_HD(:,:,run_idx{i}),rfp_HD(:,:,run_idx{i}),settings.fs,[-5 5],masks.*vessel_mask,1);
                    iR_metadata.xcorr.lag = allen.lag;
                    iR_metadata.xcorr.allen.gfp_rfp = allen.xcorr;
                    iR_metadata.xcorr.overall.gfp_rfp = overall.xcorr;
                    img_name = sprintf('xcorr_gfprfp_Run%02i_iRun%02i',Run,i);
                    exportgraphics(f,fullfile(files.images,[img_name '.png']),'Resolution',300,'BackgroundColor',[1 1 1]);
                    savefig(f,fullfile(files.images,[img_name '.fig'])); close(f);
                end
            end
            if rfp_exist
                [allen,overall,f] = f_hemLag_dT_allen(HbT(:,:,run_idx{i}),rfp_HD(:,:,run_idx{i}),settings.fs,[-5 5],masks.*vessel_mask,1);
                iR_metadata.xcorr.lag = allen.lag;
                iR_metadata.xcorr.allen.HbT_rfp = allen.xcorr;
                iR_metadata.xcorr.overall.HbT_rfp = overall.xcorr;
                img_name = sprintf('xcorr_HbTrfp_Run%02i_iRun%02i',Run,i);
                exportgraphics(f,fullfile(files.images,[img_name '.png']),'Resolution',300,'BackgroundColor',[1 1 1]);
                savefig(f,fullfile(files.images,[img_name '.fig'])); close(f);
            end
        end
        if i == 1
            metadata = iR_metadata;
        else
            metadata(i) = iR_metadata;
        end
    end

    % save summary data
    save(fullfile(files.save_dir,sprintf('Run%02i',Run),'processing_metadata.mat'),'metadata');
end

set(0,'DefaultFigureVisible','on');