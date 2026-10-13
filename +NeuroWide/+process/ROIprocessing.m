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
% INPUTS: ROIprocessing(Mouse,Date,Runs,_)
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

function ROIprocessing(Mouse,Date,Runs,varargin)
%% parse inputs
p = inputParser;
addParameter(p,'behCam',true);
addParameter(p,'newROIs',false);
addParameter(p,'extraProcessing',true);
addParameter(p,'GRAB','none');
addParameter(p,'load_dir','/projectnb/devorlab/bcraus/HRF/1P');
addParameter(p,'save_dir','/projectnb/devorlab/bcraus/HRF/1P');
addParameter(p,'plotVisible',false);
addParameter(p,'allenRef',0);
addParameter(p,'num_cores',4);

parse(p,varargin{:});
%%
% delete(gcp('nocreate'));
% parpool('local',p.Results.num_cores);

%% organize files

files = struct;
files.load_dir = p.Results.load_dir;
files.save_dir = p.Results.save_dir;
files.load_dir = fullfile(files.load_dir,Date,Mouse);
files.save_dir = fullfile(files.save_dir,Date,Mouse,'DataAnalysis');
files.dataIn = fullfile(files.load_dir,'dataIn');
files.triggers = fullfile(files.load_dir,'Triggers');
files.images = fullfile(files.save_dir,'Images');
files.processed = fullfile(files.load_dir,'processed');
files.onephoton = fullfile(files.load_dir,'onephoton');
files.camera = fullfile(files.load_dir,'camera');

[~,~,~,] = mkdir(files.save_dir);
[~,~,~] = mkdir(files.images);

% check dataIn filetype and load
if exist([files.dataIn, '.mat'], 'file')
    files.dataIn = [files.dataIn, '.mat'];
    dataIn = load(files.dataIn);
    dataIn = dataIn.dataIn;
elseif exist([files.dataIn, '.json'], 'file')
    files.dataIn = [files.dataIn, '.json'];
    dataIn = jsondecode(fileread(files.dataIn));

    % correct dataIn.template dimensions
    for i = 1:numel(dataIn)
        dataIn(i).template = permute(dataIn(i).template,[2,3,1]);
    end
end

%% check if data is stored as .dat or .bin

list_raw_runs = dir(files.onephoton);
list_raw_runs = {list_raw_runs(3:end).name};

if exist(fullfile(files.onephoton,list_raw_runs{1},'data.bin'),'file')
    data_bin = true;
else
    data_bin = false;
end

%% find runs
if data_bin
    list_runs = dir(files.onephoton);
    list_runs = {list_runs(3:end).name};
    list_runs = cellfun(@(f) str2double(f(4:5)),list_runs);
else
    list_runs = dir(files.processed);
    list_runs = {list_runs(3:end).name};
    list_runs = cellfun(@(f) str2double(f(4:7)),list_runs);
end

if ~isempty(Runs)
    list_runs = intersect(list_runs,Runs);
end

cam_folder = fullfile(files.camera,sprintf('Run%02i',list_runs(1)));
if exist(cam_folder, 'dir')
    beh_mp4 = false;
elseif exist([cam_folder,'.tif'], 'file')
    beh_mp4 = true;
end

%% create brain and behavior masks

ROI_exist = exist(fullfile(files.save_dir,'ROIs.mat'),'file');

dataIn_idx = [dataIn.runnum]==list_runs(1);

if p.Results.newROIs || ~ROI_exist
    [brain_mask,hem] = NeuroWide.process.createROIs(dataIn(dataIn_idx).template(:,:,1));
    
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
        vessel_mask = NeuroWide.process.maskVessel(dataIn(dataIn_idx).template(:,:,idx_525),sensitivity,brain_mask);
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
        if exist(cam_folder, 'dir')

            cam_images = dir(cam_folder);cam_images(1:2) = [];
            cam_images = fullfile({cam_images.folder},{cam_images.name});

            imgPath = cam_images{round((numel(cam_images))/2)};

            t = Tiff(imgPath,'r');
            imageData = im2uint8(read(t));

            behaviorROIs = NeuroWide.process.drawBehaviorROIs(imageData);
            clear cam_folder cam_images imgPath t imageData;

        elseif exist([cam_folder,'.tif'], 'file')
            cam_folder = [cam_folder,'.tif'];
            
            n_frames = numel(imfinfo(cam_folder));
            imageData = imread(cam_folder, round(n_frames / 2));
            
            behaviorROIs = NeuroWide.process.drawBehaviorROIs(imageData);
            clear cam_folder imageData
        end
        
        save(fullfile(files.save_dir,'ROIs.mat'),'behaviorROIs','-append');        
    end
else
    load(fullfile(files.save_dir,'ROIs.mat'));
end

%% create allen masks

reg_filename = fullfile(Mouse, 'DataAnalysis/ROIs.mat');
reg_files = NeuroWide.utilities.find_reg(reg_filename, 'parcellation', p.Results.load_dir, 1);

allen_exist = who('-file', fullfile(files.save_dir, 'ROIs.mat'));
allen_exist = ismember('parcellation', allen_exist);

if p.Results.newROIs || ~allen_exist
    LEDs = [dataIn(dataIn_idx).led(:).type];

    if ~isempty(reg_files)
        allen_path = reg_files{1};
        dataIn_moving = jsondecode(fileread(fullfile(fileparts(fileparts(allen_path)), 'dataIn.json')));
        
        LEDs_moving = [dataIn_moving(1).led(:).type];
        moving = squeeze(dataIn_moving(1).template(find(LEDs_moving==525),:,:));
        fixed = squeeze(dataIn(1).template(:,:,find(LEDs==525)));

        tform = NeuroWide.process.brain_register(moving, fixed);
        
        parcellation = load(allen_path,'parcellation').parcellation;
        reg_parcells = imwarp(cat(3, parcellation.Masks(:,:,:,1), parcellation.Masks(:,:,:,2)), tform, OutputView=imref2d(size(fixed)));
        parcellation.Masks = cat(4,reg_parcells(:,:,1:12), reg_parcells(:,:,13:24));
        
        rotation = dataIn(1).rotation;

    else
    
        if p.Results.allenRef
            rfp_idx = p.Results.allenRef;
        else
            rfp_idx = find(LEDs==565);
            
            if isempty(rfp_idx)
                rfp_idx = 1;
            end
        end
    
        if data_bin
            data = NeuroWide.io.data_1P(fullfile(files.onephoton,sprintf('Run%02i',list_runs(1))));
            rotation = data.meta.rotation;
            data = std(data.raw_data(:,:,:,rfp_idx),0,3);
    
            parcellation = NeuroWide.allenAtlas.allenAtlas(data,brain_mask,hem);
        else
    
            importSettings.mp = 1;
            importSettings.cores= 4;
            importSettings.doNorm = 0;
            importSettings.doCat=0;
            importSettings.nImport = 0;
            importSettings.nChannels = numel(LEDs);
            
            data = NeuroWide.io.andorDATImporter(fullfile(files.onephoton,sprintf('Run%02i',list_runs(1))),importSettings);
            data = std(data(:,:,:,rfp_idx),0,3);
        
            rotation = dataIn(1).rotation;
            
            parcellation = NeuroWide.allenAtlas.allenAtlas(imrotate(data,rotation),brain_mask,hem);
        end
    end
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
    LEDs = [dataIn(dataIn_idx).led(:).type];
    rfp_exist = sum(LEDs==565);
    gfp_exist = sum(LEDs==470);
    
    if data_bin
        files.bin = fullfile(files.onephoton,sprintf('Run%02i',Run));
        data = NeuroWide.io.data_1P(files.bin);

        masks = sum(parcellation.Masks,4);
        masks_all = cat(3, parcellation.Masks(:,:,:,1), parcellation.Masks(:,:,:,2));

        Signals = struct;
        Signals_all = struct;

        if gfp_exist
            Signals.gfp_HD = NeuroWide.process.parcellate(data.gfp_HD,masks.*vessel_mask);
            Signals.gfp = NeuroWide.process.parcellate(data.gfp,masks.*vessel_mask);
            Signals_all.gfp_HD = NeuroWide.process.parcellate(data.gfp_HD,masks_all.*vessel_mask);
            Signals_all.gfp = NeuroWide.process.parcellate(data.gfp,masks_all.*vessel_mask);
        end
        if rfp_exist
            Signals.rfp_HD = NeuroWide.process.parcellate(data.rfp_HD,masks.*vessel_mask);
            Signals.rfp = NeuroWide.process.parcellate(data.rfp,masks.*vessel_mask);
            Signals_all.rfp_HD = NeuroWide.process.parcellate(data.rfp_HD,masks_all.*vessel_mask);
            Signals_all.rfp = NeuroWide.process.parcellate(data.rfp,masks_all.*vessel_mask);
        end
    else
        files.h5 = fullfile(files.processed,sprintf('run%04i.h5',Run));
        data = struct;

        masks = sum(parcellation.Masks,4);
    
        Signals = struct;
        if gfp_exist
            data.gfp_HD = 100*h5read(files.h5,'/gfp/normHD');
            data.gfp = 100*h5read(files.h5,'/gfp/norm');
            Signals.gfp_HD = NeuroWide.process.parcellate(gfp_HD,masks.*vessel_mask);
            Signals.gfp = NeuroWide.process.parcellate(gfp,masks.*vessel_mask);
            Signals_all.gfp_HD = NeuroWide.process.parcellate(gfp_HD,masks_all.*vessel_mask);
            Signals_all.gfp = NeuroWide.process.parcellate(gfp,masks_all.*vessel_mask);
        end
        if rfp_exist
            data.rfp_HD = 100*h5read(files.h5,'/rfp/normHD');
            data.rfp = 100*h5read(files.h5,'/rfp/norm');
            Signals.rfp_HD = NeuroWide.process.parcellate(rfp_HD,masks.*vessel_mask);
            Signals.rfp = NeuroWide.process.parcellate(rfp,masks.*vessel_mask);
            Signals_all.rfp_HD = NeuroWide.process.parcellate(rfp_HD,masks_all.*vessel_mask);
            Signals_all.rfp = NeuroWide.process.parcellate(rfp,masks_all.*vessel_mask);
        end
        data.HbO = 1e6*h5read(files.h5,'/hemodynamics/HbO');
        data.HbR = 1e6*h5read(files.h5,'/hemodynamics/Hb');
        data.HbT = data.HbO+data.HbR;
    end

    Signals.HbO = NeuroWide.process.parcellate(data.HbO,masks);
    Signals.HbR = NeuroWide.process.parcellate(data.HbR,masks);
    Signals.HbT = NeuroWide.process.parcellate(data.HbT,masks);
    Signals_all.HbO = NeuroWide.process.parcellate(data.HbO,masks_all);
    Signals_all.HbR = NeuroWide.process.parcellate(data.HbR,masks_all);
    Signals_all.HbT = NeuroWide.process.parcellate(data.HbT,masks_all);

    %% process behavior

    Signals.Acc = NeuroWide.process.alignAcc(digitalInput.SonaTrigger,digitalInput.Accelerometer);
    sampling_win = floor(settings.DAQFrequency/settings.fs);
    Signals.Acc = movstd(Signals.Acc,sampling_win,0);
    Signals.Acc = Signals.Acc(round(sampling_win/2):sampling_win:end);

    Signals.t = 1/settings.fs:1/settings.fs:settings.totalTime;
    
    if p.Results.behCam
    % behavior extraction
        if beh_mp4
            cam_folder = fullfile(files.camera,sprintf('Run%02i',Run));
            cam_folder = [cam_folder,'.tif'];
            cam_images = double(NeuroWide.io.loadTifStack(cam_folder));
        else
            cam_folder = fullfile(files.camera,sprintf('Run%02i',Run));
            cam_images = dir([cam_folder,'/*.tiff']);
            cam_images = fullfile({cam_images.folder},{cam_images.name});
    
            frame_num = regexp(cam_images,'\d+','match');
            frame_num = cellfun(@(x) str2double(x{end}), frame_num);
            
            [~,order] = sort(frame_num);
            
            cam_images = cam_images(order);
        end
        % pupil
        Signals.pupil = NeuroWide.process.processPupil(cam_images,behaviorROIs);
    
        % whisker
        [Signals.whisker_long,Signals.whisker_pad] = NeuroWide.process.processWhisking(cam_images,behaviorROIs);
    
        beh_padding = ones(1,(settings.nCycles-size(Signals.pupil,2)));
        Signals.pupil = [Signals.pupil Signals.pupil(end)*beh_padding];
        Signals.whisker_long = [Signals.whisker_long Signals.whisker_long(end)*beh_padding];
        Signals.whisker_pad = [Signals.whisker_pad Signals.whisker_pad(end)*beh_padding];
    end

    % save Signals
    save(fullfile(files.save_dir,sprintf('Run%02i',Run),'Signals.mat'),'Signals','Signals_all');
    
    %% create images
    % Overview plots
    
    if ~p.Results.plotVisible
        set(0,'DefaultFigureVisible','off');
    end

    f = NeuroWide.plot.overview(Signals,p.Results.GRAB);
    exportgraphics(f,fullfile(files.images,sprintf('plotOverview_Run%02i.png',Run)),'Resolution',300,'BackgroundColor',[1 1 1]);
    savefig(f,fullfile(files.images,sprintf('plotOverview_Run%02i.fig',Run))); close(f);
        
    % Correlation images
    metadata = struct;
    if gfp_exist
        [metadata.rMaps.gfpHbT,f] = NeuroWide.connectivity.corr(brain_mask.*data.gfp_HD,data.HbT,3,plot=true);clim([-1 1]);axis image off;
        title('GFP vs. HbT');
        img_name = sprintf('gfpHbTCorr_Run%02i',Run);
        exportgraphics(f,fullfile(files.images,[img_name '.png']),'Resolution',300,'BackgroundColor',[1 1 1]);
        savefig(f,fullfile(files.images,[img_name '.fig'])); close(f);
    end

    if rfp_exist
        if gfp_exist
            [metadata.rMaps.rfpgfp,f] = NeuroWide.connectivity.corr(brain_mask.*data.gfp_HD,data.rfp_HD,3,plot=true);clim([-1 1]);axis image off;
            title('GFP vs. RFP');
            img_name = sprintf('gfprfpCorr_Run%02i',Run);
            exportgraphics(f,fullfile(files.images,[img_name '.png']),'Resolution',300,'BackgroundColor',[1 1 1]);
            savefig(f,fullfile(files.images,[img_name '.fig'])); close(f);
        end
        [metadata.rMaps.rfpHbT,f] = NeuroWide.connectivity.corr(brain_mask.*data.rfp_HD,data.HbT,3,plot=true);clim([-1 1]);axis image off;
        title('RFP vs. HbT');
        img_name = sprintf('rfpHbTCorr_Run%02i',Run);
        exportgraphics(f,fullfile(files.images,[img_name '.png']),'Resolution',300,'BackgroundColor',[1 1 1]);
        savefig(f,fullfile(files.images,[img_name '.fig'])); close(f);
    end

    if p.Results.extraProcessing        
        % calculate coherence
        if gfp_exist
            [metadata.coh.C.gfp_HbT, metadata.coh.phi.gfp_HbT, metadata.coh.fr] = NeuroWide.spectra.coherence(Signals.gfp_HD, Signals.HbT, settings.fs);
            
            f = plotCoherence(metadata.coh.C.gfp_HbT, metadata.coh.phi.gfp_HbT, metadata.coh.fr, parcellation.labelsLR, [0, 0.5]);
            title('GFP, HbT coherence');
            img_name = sprintf('coh_gfpHbT_Run%02i',Run);
            exportgraphics(f,fullfile(files.images,[img_name '.png']),'Resolution',300,'BackgroundColor',[1 1 1]);
            savefig(f,fullfile(files.images,[img_name '.fig'])); close(f);
            if rfp_exist
                [metadata.coh.C.rfp_gfp, metadata.coh.phi.rfp_gfp] = NeuroWide.spectra.coherence(Signals.rfp_HD, Signals.gfp_HD, settings.fs);
                
                f = plotCoherence(metadata.coh.C.rfp_gfp, metadata.coh.phi.rfp_gfp, metadata.coh.fr, parcellation.labelsLR, [0, 5]);
                title('RFP, GFP coherence');
                img_name = sprintf('coh_rfpgfp_Run%02i',Run);
                exportgraphics(f,fullfile(files.images,[img_name '.png']),'Resolution',300,'BackgroundColor',[1 1 1]);
                savefig(f,fullfile(files.images,[img_name '.fig'])); close(f);
            end
        end
        if rfp_exist
            [metadata.coh.C.rfp_HbT, metadata.coh.phi.rfp_HbT] = NeuroWide.spectra.coherence(Signals.rfp_HD, Signals.HbT, settings.fs);
            
            f = plotCoherence(metadata.coh.C.rfp_HbT, metadata.coh.phi.rfp_HbT, metadata.coh.fr, parcellation.labelsLR, [0, 0.5]);
            title('RFP, HbT coherence');
            img_name = sprintf('coh_rfpHbT_Run%02i',Run);
            exportgraphics(f,fullfile(files.images,[img_name '.png']),'Resolution',300,'BackgroundColor',[1 1 1]);
            savefig(f,fullfile(files.images,[img_name '.fig'])); close(f);
        end

        % calculate spectra
        if gfp_exist
            [metadata.spc.gfp, metadata.spc.fr] = NeuroWide.spectra.spectra(Signals.gfp_HD, settings.fs);
            
            f = plotSpectra(metadata.spc.gfp, metadata.spc.fr, parcellation.labelsLR, [0, 5]);
            title('GFP spectra');
            img_name = sprintf('spectra_gfp_Run%02i',Run);
            exportgraphics(f,fullfile(files.images,[img_name '.png']),'Resolution',300,'BackgroundColor',[1 1 1]);
            savefig(f,fullfile(files.images,[img_name '.fig'])); close(f);
        end
        if rfp_exist
            [metadata.spc.rfp] = NeuroWide.spectra.spectra(Signals.rfp_HD, settings.fs);
            
            f = plotSpectra(metadata.spc.rfp, metadata.spc.fr, parcellation.labelsLR, [0, 5]);
            title('RFP spectra');
            img_name = sprintf('spectra_rfp_Run%02i',Run);
            exportgraphics(f,fullfile(files.images,[img_name '.png']),'Resolution',300,'BackgroundColor',[1 1 1]);
            savefig(f,fullfile(files.images,[img_name '.fig'])); close(f);
        end
        [metadata.spc.HbT] = NeuroWide.spectra.spectra(Signals.HbT, settings.fs);
        
        f = plotSpectra(metadata.spc.HbT, metadata.spc.fr, parcellation.labelsLR, [0, 5]);
        title('HbT spectra');
        img_name = sprintf('spectra_HbT_Run%02i',Run);
        exportgraphics(f,fullfile(files.images,[img_name '.png']),'Resolution',300,'BackgroundColor',[1 1 1]);
        savefig(f,fullfile(files.images,[img_name '.fig'])); close(f);

        % calculate xcorr
        if gfp_exist
            metadata.xc.gfp_HbT = NeuroWide.IRF.xcorr(Signals.gfp_HD, Signals.HbT, 5 * settings.fs);
            metadata.xcorr.lag = -5:1/settings.fs:5;
            
            f = plotXC(metadata.xc.gfp_HbT, metadata.xcorr.lag, parcellation.labelsLR, [-5, 5]);
            title('GFP, HbT cross-correlation');
            img_name = sprintf('xcorr_HbTgfp_Run%02i',Run);
            exportgraphics(f,fullfile(files.images,[img_name '.png']),'Resolution',300,'BackgroundColor',[1 1 1]);
            savefig(f,fullfile(files.images,[img_name '.fig'])); close(f);
            if rfp_exist
                metadata.xc.rfp_gfp = NeuroWide.IRF.xcorr(Signals.rfp_HD, Signals.gfp_HD, 5 * settings.fs);
                
                f = plotXC(metadata.xc.rfp_gfp, metadata.xcorr.lag, parcellation.labelsLR, [-5, 5]);
                title('RFP, GFP cross-correlation');
                img_name = sprintf('xcorr_gfprfp_Run%02i',Run);
                exportgraphics(f,fullfile(files.images,[img_name '.png']),'Resolution',300,'BackgroundColor',[1 1 1]);
                savefig(f,fullfile(files.images,[img_name '.fig'])); close(f);
            end
        end
        if rfp_exist
            metadata.xc.rfp_HbT = NeuroWide.IRF.xcorr(Signals.rfp_HD, Signals.HbT, 5 * settings.fs);
            
            f = plotXC(metadata.xc.rfp_HbT, metadata.xcorr.lag, parcellation.labelsLR, [-5, 5]);
            title('RFP, HbT cross-correlation');
            img_name = sprintf('xcorr_HbTrfp_Run%02i',Run);
            exportgraphics(f,fullfile(files.images,[img_name '.png']),'Resolution',300,'BackgroundColor',[1 1 1]);
            savefig(f,fullfile(files.images,[img_name '.fig'])); close(f);
        end
    end

    % STIMULATION
    if settings.StimulusRepetitions
        LED = find(diff(digitalInput.LED470Signal)==1);
        LEDstim= find(diff(digitalInput.StimulusTrigger)==1);
        LEDstim = LEDstim(1:settings.StimulusDurationLong*settings.StimulusFrequency:end);
        [~,idx] = min(abs(LEDstim' - LED));

        rfp_stim = NeuroWide.process.trialAvg(Signals.rfp_HD, idx, 50, 100, dim=1, rem_baseline=1);
        gfp_stim = NeuroWide.process.trialAvg(Signals.gfp_HD, idx, 50, 100, dim=1, rem_baseline=1);
        HbT_stim = NeuroWide.process.trialAvg(Signals.HbT, idx, 50, 100, dim=1, rem_baseline=1);
        HbO_stim = NeuroWide.process.trialAvg(Signals.HbO, idx, 50, 100, dim=1, rem_baseline=1);
        HbR_stim = NeuroWide.process.trialAvg(Signals.HbR, idx, 50, 100, dim=1, rem_baseline=1);
        rfp_map = NeuroWide.process.trialAvg(data.rfp_HD, idx, 50, 100, dim=3, rem_baseline=1);
        gfp_map = NeuroWide.process.trialAvg(data.gfp_HD, idx, 50, 100, dim=3, rem_baseline=1);
        HbT_map = NeuroWide.process.trialAvg(data.HbT, idx, 50, 100, dim=3, rem_baseline=1);

        if settings.doVIS
            stim_type = 'Visual';
            reg = 12;
        else
            stim_type = 'Air Puff';
            reg = 3;
        end
        
        f = figure(Position=[100,100,1000,700]);
        tl = tiledlayout(3,3,TileSpacing='compact',Padding='compact');
        title(tl, sprintf([stim_type, ' Stimulation. %01i Hz for %01i s'], settings.StimulusFrequency, settings.StimulusDurationLong), FontSize=16);
        
        nexttile(1);
        plot(-5:0.1:10,rfp_stim(:,reg) / std(rfp_stim(:,reg), 0),color=NeuroWide.colors.r);
        xregion(0,2);
        box off;
        AX = gca;
        AX.XAxis.Visible = 'off';
        ylabel('rfp (Z)');
        
        nexttile(4);
        plot(-5:0.1:10,gfp_stim(:,reg) / std(gfp_stim(:,reg), 0),color=NeuroWide.colors.g);
        xregion(0,2);
        box off;
        AX = gca;
        AX.XAxis.Visible = 'off';
        ylabel('gfp (Z)');
        
        nexttile(7);
        hold on;
        plot(-5:0.1:10,HbT_stim(:,reg) / std(HbT_stim(:,reg), 0),color=[0,0,0]);
        plot(-5:0.1:10,HbO_stim(:,reg) / std(HbO_stim(:,reg), 0),color=NeuroWide.colors.r);
        plot(-5:0.1:10,HbR_stim(:,reg) / std(HbR_stim(:,reg), 0),color=NeuroWide.colors.b);
        xregion(0,2);
        xlabel('Time since stimulus (s)');
        box off;
        ylabel('Hb (Z)');
        
        rfp_ratio = (max(rfp_map(:,:,51:end), [], 3) - mean(rfp_map(:,:,1:50),3)) ./ std(rfp_map(:,:,1:50),0,3);
        gfp_ratio = (max(gfp_map(:,:,51:end), [], 3) - mean(gfp_map(:,:,1:50),3)) ./ std(gfp_map(:,:,1:50),0,3);
        HbT_ratio = (max(HbT_map(:,:,51:end), [], 3) - mean(HbT_map(:,:,1:50),3)) ./ std(HbT_map(:,:,1:50),0,3);
        
        single_cmp = cat(2, linspace(0,1,256)', zeros(256,1), zeros(256,1));
        AX = nexttile(2);
        imagesc(rfp_ratio, AlphaData=brain_mask);
        cb = colorbar;
        cb.Label.String = 'Z';
        clim([0, prctile(rfp_ratio(:), 99)]);
        axis image off;
        colormap(AX, single_cmp);
        
        AX = nexttile(5);
        imagesc(gfp_ratio, AlphaData=brain_mask);
        cb = colorbar;
        cb.Label.String = 'Z';
        clim([0, prctile(gfp_ratio(:), 99)]);
        axis image off;
        colormap(AX, single_cmp(:,[2,1,3]));
        
        AX = nexttile(8);
        imagesc(HbT_ratio, AlphaData=brain_mask);
        cb = colorbar;
        cb.Label.String = 'Z';
        clim([0, prctile(HbT_ratio(:), 99)]);
        axis image off;
        colormap(AX, single_cmp(:,[2,3,1]));
        
        rgb_img = cat(3, rfp_ratio / prctile(rfp_ratio(:), 99), ...
                        gfp_ratio / prctile(gfp_ratio(:), 99), ...
                        HbT_ratio / prctile(HbT_ratio(:), 99));
        
        nexttile(3);
        h = imshow(rgb_img);
        h.AlphaData = brain_mask;
        
        nexttile(6);
        XC_gfp_rfp = NeuroWide.IRF.xcorr(gfp_stim, rfp_stim, settings.fs*5);
        XC_gfp_HbT = NeuroWide.IRF.xcorr(gfp_stim, HbT_stim, settings.fs*5);
        XC_rfp_HbT = NeuroWide.IRF.xcorr(rfp_stim, HbT_stim, settings.fs*5);
        hold on;
        plot(-5:1/settings.fs:5,XC_gfp_rfp(:,reg), color=NeuroWide.colors.y);
        plot(-5:1/settings.fs:5,XC_gfp_HbT(:,reg), color=NeuroWide.colors.c);
        plot(-5:1/settings.fs:5,XC_rfp_HbT(:,reg), color=NeuroWide.colors.p);
        xline(0);
        ylim([-1, 1]);
        xlim([-5, 5]);
        xlabel('Delay (s)');
        ylabel('r');
        title('Cross-correlation');
        
        mask = brain_mask;
        mask(isnan(mask)) = 0;
        mask = logical(mask);
        spatial_r = corrcoef([rfp_ratio(mask),gfp_ratio(mask),HbT_ratio(mask)]);
        nexttile(9);
        bar([spatial_r(2,1), spatial_r(3,1), spatial_r(2,3)]);
        box off;
        ylabel('r');
        xticklabels({'rfp-gfp','rfp-HbT','gfp-HbT'});
        title('Spatial Correlation');
        
        img_name = sprintf('stimulation_Run%02i', Run);
        exportgraphics(f,fullfile(files.images,[img_name '.png']),'Resolution',300,'BackgroundColor',[1 1 1]);
        savefig(f,fullfile(files.images,[img_name '.fig'])); close(f);

        % store metadata
        
        metadata.stim.rfp_ratio = rfp_ratio;
        metadata.stim.gfp_ratio = gfp_ratio;
        metadata.stim.HbT_ratio = HbT_ratio;
        
        metadata.stim.t_bef = 50;
        metadata.stim.t_aft = 150;
        
        metadata.stim.stim_idx = idx;
    end
    

    % save summary data
    save(fullfile(files.save_dir,sprintf('Run%02i',Run),'processing_metadata.mat'),'metadata');
end

set(0,'DefaultFigureVisible','on');
%%

%%%%%%%%%%%%%%%% plotting functions %%%%%%%%%%%%%%%%

function f = plotXC(R, t, labels, x_lim)

f = figure(Position=[0,0,800,800]);
tiledlayout(2, 1);

ax(1) = nexttile();
imagesc(R', XData=t);
colormap(NeuroWide.cmp.bbr);
xlim(x_lim);
clim([-1, 1]);
c = colorbar(ax(1));
c.Label.String = 'r';
set(ax(1), YTick=1:12, YTickLabels=labels, XTickLabels=[]);
ax(1).TickDir = 'out';
box off;

ax(2) = nexttile;
plot(t, mean(R, 2));
ylabel('r');
xlim(x_lim);
ylim([-1, 1]);
box off;
xlabel('Lag (s)');

end

function f = plotSpectra(S, fr, labels, x_lim)

f = figure(Position=[0,0,800,800]);
tiledlayout(2, 1);

ax(1) = nexttile();
imagesc(log10(S'), XData=fr);
colormap(NeuroWide.cmp.inferno);
xlim(x_lim);
c = colorbar(ax(1));
c.Label.String = 'log_1_0(Power)';
set(ax(1), YTick=1:12, YTickLabels=labels, XTickLabels=[]);
ax(1).TickDir = 'out';
box off;

ax(2) = nexttile;
plot(fr, log10(mean(S, 2)));
ylabel('log_1_0(Power)');
xlim(x_lim);
box off;
xlabel('F (Hz)');

end

function f = plotCoherence(C, phi, fr, labels, x_lim)

f = figure(Position=[0,0,800,1000]);
tiledlayout(3, 1);

ax(1) = nexttile();
imagesc(C', XData=fr);
colormap(NeuroWide.cmp.inferno);
clim([0, 1]);
xlim(x_lim);
c = colorbar(ax(1), Ticks=[0, 1]);
c.Label.String = 'coherence';
set(ax(1), YTick=1:12, YTickLabels=labels, XTickLabels=[]);
ax(1).TickDir = 'out';
box off;

ax(2) = nexttile();
imagesc(phi', XData=fr);
colormap(ax(2), NeuroWide.cmp.bbr);
clim(pi * [-1, 1]);
xlim(x_lim);
c = colorbar(ax(2), Ticks=[-pi, pi], TickLabels={'-\pi', '\pi'});
c.Label.String = 'radians';
set(ax(2), YTick=1:12, YTickLabels=labels, XTickLabels=[]);
ax(2).TickDir = 'out';
box off;

ax(3) = nexttile;
plot(fr, mean(C, 2));
ylim([0, 1]);
ylabel('Coherence');
yticks([0, 1]);

yyaxis right;
plot(fr, mean(phi, 2));
ylim(pi * [-1, 1]);
xlim(x_lim);
box off;
xlabel('F (Hz)');
ylabel('Phase (radians)');
yticks(pi * [-1, 1]);
yticklabels({'-\pi', '\pi'});
end

%%
end