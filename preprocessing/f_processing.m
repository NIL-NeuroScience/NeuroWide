%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%                   Single photon data preprocessing
% 
% Initial processing script for widefield fluorescence and intrinsic
% imaging data. Retrieves and saves imaging metadata. Estimates dF/F for
% fluorescent imaging data, estimates changes in hemoglobin oxygenation and
% concentration, and applies hemodynamic artifact correction algorithms to
% both fluorescent channels.
% 
% Authors: Martin Thunemann, Patrick Doran, Brad Rauscher
% 
% INPUTS: f_ROIprocessing(Mouse,Date,Runs,_)
%   Mouse: mouse ID
%   Date: session date
% 
% EXTRA PARAMETERS:
%   runs: list of runs to analyze. Leave empty to analyze all available
%       runs (default = [])
%   smooth: 2D smoothing kernel. Set to 0 to not smooth (default = 2)
%   compression: compression value (default = 0)
%   rotation: rotation of videos in degrees (default = '0')
%   load_dir: directory to load from (default = 'bcraus/HRF/1P')
%   save_dir: directory to save to (default = 'bcraus/HRF/1P')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f_processing(Mouse,Date,varargin)
%% parse inputs
p = inputParser;
addParameter(p,'runs',[]);
addParameter(p,'smooth',2);
addParameter(p,'compression',0);
addParameter(p,'rotation',0);
addParameter(p,'load_dir','bcraus/HRF/1P');
addParameter(p,'save_dir','bcraus/HRF/1P');

parse(p,varargin{:});

%% Define settings for data processing

commonSettings.Runs = p.Results.runs;

commonSettings.hasPylon = false;        % True when there are pylon recordings
commonSettings.compression = p.Results.compression;
commonSettings.chunkSize = [2^4 2^4 2^4]; % optimized for L2 cache of mesoscope computer
commonSettings.smoothVal = p.Results.smooth;
commonSettings.rotation = p.Results.rotation;

%% Select root folder
rootFolderList={fullfile('/projectnb/devorlab',p.Results.load_dir,Date,Mouse)};
%%
for rootFolderCounter=1:size(rootFolderList,1)
    %%
    clearvars dataIn commonFolders
    commonFolders.root=rootFolderList{rootFolderCounter};
    if ~exist(commonFolders.root,'dir')
        error([commonFolders.root ' does not exist']);
    end

    %% Reads root folder and returns folder structure
    commonFolders.folderIn=f_returnFolderStructure(commonFolders.root,4);
    commonFolders.folderList=f_findFoldersMaster(commonFolders.folderIn,{'onephoton','daq','behaviorCam','processed','ephys'});
    
    %% Lists all runs from solis input with tif files and sorts
    clearvars tmp*
    tmpSolis=find(strcmp({commonFolders.folderIn.type},'1photon'));
    if isempty(tmpSolis);return;end
    if size(tmpSolis,2)==1
        for folderCounter=1:size(commonFolders.folderIn(tmpSolis).folders,2)
            clearvars tmp1 tmp2 tmp3 tmpNum
            commonFolders.solis(folderCounter).runnum=commonFolders.folderIn(tmpSolis).folders(folderCounter).runnum;
            commonFolders.solis(folderCounter).name=commonFolders.folderIn(tmpSolis).folders(folderCounter).name;
            if isfield(commonFolders.folderIn(tmpSolis).folders(folderCounter),'folders') && ~isempty(commonFolders.folderIn(tmpSolis).folders(folderCounter).folders)
                tmpSifx=find(strcmp(commonFolders.folderIn(tmpSolis).folders(folderCounter).files(:,2),'sifx'));
                tmpDat=find(strcmp(commonFolders.folderIn(tmpSolis).folders(folderCounter).files(:,2),'dat'));
                tmpIni=find(strcmp(commonFolders.folderIn(tmpSolis).folders(folderCounter).files(:,2),'ini'));
                if size(tmpSifx,2)==1 && strcmp(commonFolders.folderIn(tmpSolis).folders(folderCounter).files{tmpSifx,1},'Spooled files.sifx')
                    tmp0(1)=true;
                end
                if size(tmpIni,2)==1 && strcmp(commonFolders.folderIn(tmpSolis).folders(folderCounter).files{tmpIni,1},'acquisitionmetadata.ini')
                    tmp0(2)=true;
                end
            end
            if all(tmp0)
                commonFolders.solis(folderCounter).isSifx=true;
            else
                commonFolders.solis(folderCounter).isSifx=false;
            end
            commonFolders.solis(folderCounter).folder=[commonFolders.folderIn(tmpSolis).folders(folderCounter).parentFolder,filesep,commonFolders.folderIn(tmpSolis).folders(folderCounter).name];
            commonFolders.solis(folderCounter).isSona=false;
            commonFolders.solis(folderCounter).isZyla=false;
            commonFolders.solis(folderCounter).isDualCam=false;
        end
    elseif size(tmpSolis,2)>1
        for iFolder=1:size(tmpSolis,2)
            for folderCounter=1:size(commonFolders.folderIn(tmpSolis(iFolder)).folders,2)
                clearvars tmp1 tmp2 tmp3 tmpNum
                commonFolders.solis(iFolder,folderCounter).runnum=commonFolders.folderIn(tmpSolis(iFolder)).folders(folderCounter).runnum;
                commonFolders.solis(iFolder,folderCounter).name=commonFolders.folderIn(tmpSolis(iFolder)).folders(folderCounter).name;
                if isfield(commonFolders.folderIn(tmpSolis(iFolder)).folders(folderCounter),'folders') && ~isempty(commonFolders.folderIn(tmpSolis(iFolder)).folders(folderCounter).folders)
                    tmpSifx=find(strcmp(commonFolders.folderIn(tmpSolis(iFolder)).folders(folderCounter).files(:,2),'sifx'));
                    tmpDat=find(strcmp(commonFolders.folderIn(tmpSolis(iFolder)).folders(folderCounter).files(:,2),'dat'));
                    tmpIni=find(strcmp(commonFolders.folderIn(tmpSolis(iFolder)).folders(folderCounter).files(:,2),'ini'));
                    if size(tmpSifx,2)==1 && strcmp(commonFolders.folderIn(tmpSolis(iFolder)).folders(folderCounter).files{tmpSifx,1},'Spooled files.sifx')
                        tmp0(1)=true;
                    end
                    if size(tmpIni,2)==1 && strcmp(commonFolders.folderIn(tmpSolis(iFolder)).folders(folderCounter).files{tmpIni,1},'acquisitionmetadata.ini')
                        tmp0(2)=true;
                    end
                end
                if all(tmp0)
                    commonFolders.solis(iFolder,folderCounter).isSifx=true;
                else
                    commonFolders.solis(iFolder,folderCounter).isSifx=false;
                end
                commonFolders.solis(iFolder,folderCounter).folder=[commonFolders.folderIn(tmpSolis(iFolder)).folders(folderCounter).parentFolder,filesep,commonFolders.folderIn(tmpSolis(iFolder)).folders(folderCounter).name];

                if contains(commonFolders.folderIn(tmpSolis(iFolder)).folders(folderCounter).parentFolder,'Sona','IgnoreCase',true) && ~contains(commonFolders.folderIn(tmpSolis(iFolder)).folders(folderCounter).parentFolder,'Zyla','IgnoreCase',true)
                    commonFolders.solis(iFolder,folderCounter).isSona=true;
                    commonFolders.solis(iFolder,folderCounter).isZyla=false;
                    commonFolders.solis(iFolder,folderCounter).isDualCam=true;
                elseif contains(commonFolders.folderIn(tmpSolis(iFolder)).folders(folderCounter).parentFolder,'Zyla','IgnoreCase',true) && ~contains(commonFolders.folderIn(tmpSolis(iFolder)).folders(folderCounter).parentFolder,'Sona','IgnoreCase',true)
                    commonFolders.solis(iFolder,folderCounter).isZyla=true;
                    commonFolders.solis(iFolder,folderCounter).isSona=false;
                    commonFolders.solis(iFolder,folderCounter).isDualCam=true;
                else
                    commonFolders.solis(iFolder,folderCounter).isSona=false;
                    commonFolders.solis(iFolder,folderCounter).isZyla=false;
                end
            end
        end
    end

    clearvars i* tmp*
    
    %% Lists all mat files generated with DAQ/Triggermaster control code
    tmpDAQ=find(strcmp({commonFolders.folderIn.type},'daq'));
    tmpEphysIdx = find(strcmp({commonFolders.folderIn.type},'ephys'));
    if isempty(tmpDAQ);return;end
    if isfield(commonFolders.folderIn(tmpDAQ),'files') && ~isempty(commonFolders.folderIn(tmpDAQ).files)
        for iFile=1:size(commonFolders.folderIn(tmpDAQ).files,1)
            if ~isempty(commonFolders.folderIn(tmpDAQ).files(iFile,3))
                commonFolders.daq(iFile).runnum=commonFolders.folderIn(tmpDAQ).files{iFile,3};
                commonFolders.daq(iFile).folder=[commonFolders.folderIn(tmpDAQ).parentFolder,filesep,commonFolders.folderIn(tmpDAQ).name];
                commonFolders.daq(iFile).name=commonFolders.folderIn(tmpDAQ).files{iFile,1};
            else
                commonFolders.daq(iFile).runnum=NaN;
            end
            
            if ~isempty(tmpEphysIdx) && ~isempty(commonFolders.folderIn(tmpEphysIdx).folders(iFile))
                ephysFileIdx = find([commonFolders.folderIn(tmpEphysIdx).folders.runnum]==commonFolders.daq(iFile).runnum);
                commonFolders.ephys(iFile).runnum=commonFolders.folderIn(tmpDAQ).files{iFile,3};
                commonFolders.ephys(iFile).folder=fullfile(commonFolders.folderIn(tmpEphysIdx).folders(ephysFileIdx).parentFolder,commonFolders.folderIn(tmpEphysIdx).folders(ephysFileIdx).name);
                commonFolders.ephys(iFile).name=commonFolders.folderIn(tmpEphysIdx).folders(ephysFileIdx).files;
            else
                commonFolders.ephys(iFile).runnum=NaN;
            end
        end
    end

    %% Generates dataIn variable with individual runs and references DAQ/behavior recording to the runs
    tmpTable=[];
    tmpEphysIdx = strcmp({commonFolders.folderIn.type},'ephys');
    for iEntry=1:size(commonFolders.solis,1)
        for iEntry2=1:size(commonFolders.solis,2)
        if commonFolders.solis(iEntry,iEntry2).isSifx
            tmpFind=find([commonFolders.daq.runnum]==commonFolders.solis(iEntry,iEntry2).runnum);
            if ~isempty(tmpFind) && size(tmpFind,2)==1
                tmpTable(end+1).runnum=commonFolders.solis(iEntry,iEntry2).runnum;
                tmpTable(end).solis=commonFolders.solis(iEntry,iEntry2);
                tmpTable(end).isSona=commonFolders.solis(iEntry,iEntry2).isSona;
                tmpTable(end).isZyla=commonFolders.solis(iEntry,iEntry2).isZyla;
                tmpTable(end).isDualCam=commonFolders.solis(iEntry,iEntry2).isDualCam;
                tmpTable(end).daq=commonFolders.daq(tmpFind);
                if ~isempty(tmpEphysIdx)
                    tmpTable(end).ephys=commonFolders.ephys(tmpFind);
                end
            else
                tmpTable(end+1).daqID=[];
            end
            if commonSettings.hasPylon
                tmpFind=find([commonFolders.pylon.runnum]==commonFolders.solis(iEntry).runnum);
                if ~isempty(tmpFind) && size(tmpFind,2)==1
                    tmpTable(end).pylon=commonFolders.pylon(tmpFind);
                else
                    tmpTable(end).pylon=[];
                end
            end
        end
        end
    end
    dataIn=tmpTable;
    clearvars tmp* i*

    %% Loads parameters from DAQ/Triggermaster control file and performs consistency checks
    for folderCounter=1:size(dataIn,2)
        tmpDAQFile=[dataIn(folderCounter).daq.folder,filesep,dataIn(folderCounter).daq.name];
        tmpWhos=whos('-file',tmpDAQFile); %this is slow, better way of doing this?
        if any(strcmp({tmpWhos.name},'settings'))
            tmpIn=load(tmpDAQFile,'settings');
            tmpDigIn=load(tmpDAQFile,'digitalInput');

            dataIn(folderCounter).settings=tmpIn.settings;
            dataIn(folderCounter).frameNumberConistent=true;
            if commonSettings.hasPylon
                if ~isempty(dataIn(folderCounter).pylon) && dataIn(folderCounter).settings.nCycles==dataIn(folderCounter).pylon.tifN
                    dataIn(folderCounter).cycleNumberConistent=true;
                else
                    dataIn(folderCounter).cycleNumberConistent=false;
                end
            end
        else
            dataIn(folderCounter).settings=[];
        end
    end
    clearvars tmp* i*

    %% Add quality control
    for folderCounter=1:size(dataIn,2)
        if dataIn(folderCounter).frameNumberConistent
            dataIn(folderCounter).goodRun=true;
        else
            dataIn(folderCounter).goodRun=false;
        end
    end

    if ~isempty(commonSettings.Runs)
        for iRun=1:size(dataIn,2)
            if ~ismember(dataIn(iRun).runnum,commonSettings.Runs)
                dataIn(iRun).goodRun=false;
            end
        end
    end
    clearvars tmp* i*

    %% Matches frames with LED and loads some images as template
    tmpSettings.mp = 1;tmpSettings.cores= 4;tmpSettings.doNorm = 0;tmpSettings.doCat=0;
    tmpSettings.nImport = 0;
    
    for folderCounter=find([dataIn.goodRun])
        tmpSettings.frames = dataIn(folderCounter).settings.nframes;
        tmpSettings.nChannels=size(dataIn(folderCounter).settings.LEDOrder,1);
        tmpRaw=f_AndorDATImporter(dataIn(folderCounter).solis.folder,tmpSettings);
        for iLED=1:size(dataIn(folderCounter).settings.LEDOrder,1)
            dataIn(folderCounter).led(iLED).type=str2double(dataIn(folderCounter).settings.LEDOrder(iLED,:));
            dataIn(folderCounter).led(iLED).time=dataIn(folderCounter).settings.ExposureTimes(iLED);
            dataIn(folderCounter).led(iLED).power=dataIn(folderCounter).settings.LEDPower(iLED);

            dataIn(folderCounter).template(:,:,iLED)=mean(tmpRaw(:,:,:,iLED),3);
        end
        dataIn(folderCounter).template = imrotate(dataIn(folderCounter).template,commonSettings.rotation);
        dataIn(folderCounter).rotation = commonSettings.rotation;
    end
    clearvars tmp* i*

    %% Generates folder 'processed' and folder 'images' for output

    commonFolders.processed=fullfile('/projectnb/devorlab',p.Results.save_dir,Date,Mouse,'processed');
    [~,~,~] = mkdir(commonFolders.processed);

    clearvars tmp* i*

    %% Store dataIn variable with information in root folder
    save(fullfile('/projectnb/devorlab',p.Results.save_dir,Date,Mouse,'dataIn.mat'),'dataIn');
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Preprocessing is finished. Start the image processing run-by-run
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%
    for folderCounter=find([dataIn.goodRun])

        if exist(fullfile(commonFolders.processed,sprintf('run%04i.h5',dataIn(folderCounter).runnum)),'file')
            delete(fullfile(commonFolders.processed,sprintf('run%04i.h5',dataIn(folderCounter).runnum)));
        end

        %%
        clearvars settings* tmp* gfp* rfp* Hb*
        settings=dataIn(folderCounter);
        settings.hdChannel1=find([dataIn(folderCounter).led.type]==525);
        settings.hdChannel2=find([dataIn(folderCounter).led.type]==625);
        settings.gfpChannel=find([dataIn(folderCounter).led.type]==470);
        settings.rfpChannel=find([dataIn(folderCounter).led.type]==565);
        settings.hd=f_defineHemodynamicParameters(525,625);

        %% Importing the entire image dataset
        tmpSettings.mp = 1;tmpSettings.cores= 4;tmpSettings.doNorm = 0;tmpSettings.nChannels = size(dataIn(folderCounter).settings.LEDOrder,1);tmpSettings.doCat=0;tmpSettings.nImport=0;
        tmpSettings.frames = dataIn(folderCounter).settings.nframes;
        rawImageChannel=f_AndorDATImporter(dataIn(folderCounter).solis.folder,tmpSettings);
        
        if ~isempty(settings.rfpChannel)
            rfpChannel = imrotate(rawImageChannel(:,:,:,settings.rfpChannel),commonSettings.rotation);
        end
        if ~isempty(settings.gfpChannel)
            gfpChannel = imrotate(rawImageChannel(:,:,:,settings.gfpChannel),commonSettings.rotation);
        end
        if ~isempty(settings.hdChannel1)
            hdChannel1 = imrotate(rawImageChannel(:,:,:,settings.hdChannel1),commonSettings.rotation);
            hdChannel2 = imrotate(rawImageChannel(:,:,:,settings.hdChannel2),commonSettings.rotation);
        end

        clearvars tmpSettings rawImageChannel;
        
        % PROCESSING BEGINS

        %% Process reflectance images to estimate [Hb]
        if ~isempty(settings.hdChannel1)
            fprintf('\nCalculate Hb, HbO...')

            [Hb,HbO] = f_calcHb(hdChannel1,hdChannel2);
            clearvars tmp* i*
            
            HbO = process_video(HbO,commonSettings.smoothVal);
            Hb = process_video(Hb,commonSettings.smoothVal);
            
            fprintf('\nSaving Hb and HbO data (h5)...')
            save_h5(Hb,fullfile(commonFolders.processed,sprintf('run%04i.h5',dataIn(folderCounter).runnum)),'/hemodynamics/Hb',commonSettings.compression,commonSettings.chunkSize);
            save_h5(HbO,fullfile(commonFolders.processed,sprintf('run%04i.h5',dataIn(folderCounter).runnum)),'/hemodynamics/HbO',commonSettings.compression,commonSettings.chunkSize);
            fprintf('done.\n');
        end
        
        %% process gfp
        if ~isempty(settings.gfpChannel)
            fprintf('\nProcess green fluorescence...')
            
            [gfp_norm,gfp_norm_HD] = f_470HD(gfpChannel,HbO,Hb);

            gfp_norm = process_video(gfp_norm,commonSettings.smoothVal);
            gfp_norm_HD = process_video(gfp_norm_HD,commonSettings.smoothVal);

            clearvars tmp* i*
            fprintf('\nSaving green fluorescence data (h5)...')
            save_h5(gfp_norm,fullfile(commonFolders.processed,sprintf('run%04i.h5',dataIn(folderCounter).runnum)),'/gfp/norm',commonSettings.compression,commonSettings.chunkSize);
            save_h5(gfp_norm_HD,fullfile(commonFolders.processed,sprintf('run%04i.h5',dataIn(folderCounter).runnum)),'/gfp/normHD',commonSettings.compression,commonSettings.chunkSize);
            fprintf('done.\n');
        end
        
        %% process rfp
        if ~isempty(settings.rfpChannel)
            fprintf('\nProcess red fluorescence...')
            
            [rfp_norm,rfp_norm_HD] = f_565HD(rfpChannel,hdChannel1,hdChannel2);
            
            rfp_norm = process_video(rfp_norm,commonSettings.smoothVal);
            rfp_norm_HD = process_video(rfp_norm_HD,commonSettings.smoothVal);

            clearvars tmp* i
            fprintf('\nSaving red fluorescence data (h5)...');
            save_h5(rfp_norm,fullfile(commonFolders.processed,sprintf('run%04i.h5',dataIn(folderCounter).runnum)),'/rfp/norm',commonSettings.compression,commonSettings.chunkSize);
            save_h5(rfp_norm_HD,fullfile(commonFolders.processed,sprintf('run%04i.h5',dataIn(folderCounter).runnum)),'/rfp/normHD',commonSettings.compression,commonSettings.chunkSize);
            fprintf('done.\n');
        end
    end
end

end

function processed_signal = process_video(signal,smooth)
    signal(isnan(signal)) = 0;
    if smooth
        processed_signal = f_smooth2d(signal,smooth);
    end
end

function save_h5(signal,filename,path,compression,chunk_size)
    h5create(filename,path,size(signal),'Deflate',compression,'Chunksize',chunk_size,'Datatype','double');
    h5write(filename,path,signal);
end