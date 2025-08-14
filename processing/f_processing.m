%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%                   Single photon data preprocessing
% 
% Version 0.1 MTH 12/16/2021
% Version 0.2 PD  01/28/2022
% Version 0.3 PD  03/31/2022
% Version 1.0 MTH 01/28/2023
% Version 1.1 MTH 03/20/2023 processing of multiple files/deactivated mask
% Version 1.2 MTH 03/20/2023 rearranged save and clearvars to avoid error 
%                            for very large files
% Version 1.3 MTH 04/26/2023 improved handling of series w/o 4 channels
%                            turned off try catch for debugging
% Version 1.4 BCR 07/11/2023 added smoothing option and specific baseline
%                            option, added rfp HD correction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f_processing(date,mouse,rotation,allRuns,dataDir,user,badRuns)
%% test inputs
% date = '25-06-10';
% mouse = 'B6_247';
% rotation = 180;
% allRuns = 1;
% dataDir = 'HRF/1P';
% user = 'bcraus';
% nargin = 6;

%% Adds folders with functions/subfunctions to the MATLAB search path
commonSettings.matlabCodeRoot='/projectnb/devorlab/bcraus';

% addpath('/projectnb/devorlab/bcraus/AnalysisCode/1P/Single_Photon_Processing_Digital')
% addpath('/projectnb/devorlab/bcraus/AnalysisCode/1P/Single_Photon_Processing_MTH2023')
% addpath('/projectnb/devorlab/bcraus/AnalysisCode/processing_functions');

%% Define settings for data processing
%commonSettings.binFactor=1;                % binning of the images (1=off;try 2 or 4)
commonSettings.templateSize=10;             % number of frames loaded to generate template image
commonSettings.maskLED=470;                 % images with this LED are shown to define mask for image cropping
%commonSettings.tolerance = 0.6;            % Portion of region that has to be in image to extract timecourse
commonSettings.allRuns = allRuns;                 % Set to 1 to do all runs or 0 to only do certain runs
if nargin > 6
    commonSettings.badRuns = badRuns;         % Vector of runs not To Process. Only used if All runs is 0
else
    commonSettings.badRuns = [];
end

commonSettings.hasPylon = false;            % True when there are pylon recordings
commonSettings.saveHDF5 = true;         % Stores data as h5 file
commonSettings.channelList = [1,1]; % F-green, F-red, R-530, R-630; to be implemented.
commonSettings.debug.txt='/projectnb/devorlab/bcraus/fullProcessing.txt';% verbose output, to be implemented.
commonSettings.compression=0;
commonSettings.chunkSize=[2^4 2^4 2^4];          % optimized for L2 cache of mesoscope computer
commonSettings.gfprfpswitched = false;      % sometimes gfp and rfp are switched, set to true if this is the case               
commonSettings.smoothing = true;
commonSettings.smoothVal = 2;
commonSettings.rotation = rotation;

%% Select root folder
rootFolderList={fullfile('/projectnb/devorlab',user,dataDir,date,mouse)};
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

    genus = zeros(size(dataIn,2),1);
    gen = struct;
    %% Loads parameters from DAQ/Triggermaster control file and performs consistency checks
    for folderCounter=1:size(dataIn,2)
        tmpDAQFile=[dataIn(folderCounter).daq.folder,filesep,dataIn(folderCounter).daq.name];
        tmpWhos=whos('-file',tmpDAQFile); %this is slow, better way of doing this?
        if any(strcmp({tmpWhos.name},'settings'))
            tmpIn=load(tmpDAQFile,'settings');
            tmpDigIn=load(tmpDAQFile,'digitalInput');

            genus(folderCounter) = 0;%tmpIn.settings.doGENUS;
            if genus(folderCounter)
                
                gen.nChannels = size(tmpIn.settings.ExposureTimes,2);
                gen.frames = tmpIn.settings.nCycles;
                
                for nIdx = 1:gen.nChannels
                    
                    gen.rising470 = find(diff(tmpDigIn.digitalInput.(['LED' tmpIn.settings.LEDOrder(nIdx,:) 'Signal']))==1);
                    gen.falling470 = find(diff(tmpDigIn.digitalInput.(['LED' tmpIn.settings.LEDOrder(nIdx,:) 'Signal']))==-1);
                    gen.framelength470 = gen.falling470-gen.rising470;

                    for gIdx = 1:gen.frames
                        gen.r(gIdx,nIdx,folderCounter) = sum(tmpDigIn.digitalInput.StimulusTrigger(gen.rising470(gIdx):gen.falling470(gIdx)))/(gen.framelength470(gIdx)+1);
                    end
                end
            end

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

    if ~commonSettings.allRuns
        for iRun=1:size(dataIn,2)
            if ismember(dataIn(iRun).runnum,commonSettings.badRuns)
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
        tmpRaw=f_AndorDATImporter_WF(dataIn(folderCounter).solis.folder,tmpSettings);
        for iLED=1:size(dataIn(folderCounter).settings.LEDOrder,1)
            dataIn(folderCounter).led(iLED).type=str2double(dataIn(folderCounter).settings.LEDOrder(iLED,:));
            dataIn(folderCounter).led(iLED).time=dataIn(folderCounter).settings.ExposureTimes(iLED);
            dataIn(folderCounter).led(iLED).power=dataIn(folderCounter).settings.LEDPower(iLED);

            if genus(folderCounter)
                tmpRaw = double(tmpRaw);
                rMat = f_hemRegressGEN(tmpRaw(:,:,:,iLED),gen.r(:,iLED,folderCounter));
                reg = permute(gen.r(:,iLED,folderCounter),[2 3 1]);
                tmpRaw(:,:,:,iLED) = tmpRaw(:,:,:,iLED)-rMat.*reg;
                dataIn(folderCounter).rMat(:,:,iLED) = rMat;
                dataIn(folderCounter).reg = gen.r;
            end

            dataIn(folderCounter).template(:,:,iLED)=mean(tmpRaw(:,:,:,iLED),3);
        end
        dataIn(folderCounter).template = imrotate(dataIn(folderCounter).template,commonSettings.rotation);
        dataIn(folderCounter).rotation = commonSettings.rotation;
    end
    clearvars tmp* i*

    %% Generates folder 'processed' and folder 'images' for output
    tmpProcessed=find(strcmp({commonFolders.folderIn.name},'processed'));
    if isempty(tmpProcessed)
        [~,~]=mkdir([commonFolders.root,filesep,'processed']);
    end

    commonFolders.processed=fullfile('/projectnb/devorlab/bcraus',dataDir,date,mouse,'processed');
    [~,~,~] = mkdir(commonFolders.processed);
%         tmpImages=find(strcmp({commonFolders.folderIn.name},'images'));
%         if isempty(tmpImages)
%             [~,~]=mkdir([commonFolders.root,filesep,'images']);
%         end
%         commonFolders.images=[commonFolders.root,filesep,'images'];
    clearvars tmp* i*

    %% Store dataIn variable with information in root folder
    save([fullfile('/projectnb/devorlab/bcraus',dataDir,date,mouse) filesep 'dataIn.mat'],'dataIn')
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Preprocessing is finished. Start the image processing run-by-run
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%
    for folderCounter=find([dataIn.goodRun])
        if dataIn(folderCounter).isDualCam;disp('Sorry, this code does not work with dual-camera recordings!');return;end

        if commonSettings.saveHDF5
            if exist([commonFolders.processed,filesep,'run',num2str(dataIn(folderCounter).runnum,'%04.0f'),'.h5'],'file')
                delete([commonFolders.processed,filesep,'run',num2str(dataIn(folderCounter).runnum,'%04.0f'),'.h5']);
            end
        end

        %%
        clearvars settings* tmp* gfp* rfp* Hb*
        tic;
        settings=dataIn(folderCounter);
        settings.hdChannel1=find([dataIn(folderCounter).led.type]==525);
        settings.hdChannel2=find([dataIn(folderCounter).led.type]==625);
        settings.gfpChannel=find([dataIn(folderCounter).led.type]==470);
        settings.rfpChannel=find([dataIn(folderCounter).led.type]==565);
        settings.hd=f_defineHemodynamicParameters_WF(525,625);

        %% Importing the entire image dataset
%             folder = '/projectnb/devorlab/bcraus/Tests/1P/DualCam/23-06-30/crosstalk/mApple/onephotonSona/run10';
        tmpSettings.mp = 1;tmpSettings.cores= 4;tmpSettings.doNorm = 0;tmpSettings.nChannels = size(dataIn(folderCounter).settings.LEDOrder,1);tmpSettings.doCat=0;tmpSettings.nImport=0;
        tmpSettings.frames = dataIn(folderCounter).settings.nframes;
        rawImageChannel=f_AndorDATImporter_WF(dataIn(folderCounter).solis.folder,tmpSettings);
        
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
        toc;
% PROCESSING BEGINS

        %% Process reflectance images to estimate [Hb]
        if ~isempty(settings.hdChannel1)
            fprintf('\nCalculate Hb, HbO...')

            %optional: add filter and other pp

            tmpA0Hb=settings.hd.cLambda2Hb*log(mean(hdChannel2,3))-settings.hd.cLambda1Hb*log(mean(hdChannel1,3));
            tmpA0HbO=settings.hd.cLambda2HbO*log(mean(hdChannel2,3))-settings.hd.cLambda1HbO*log(mean(hdChannel1,3));
            HbO=(tmpA0HbO+settings.hd.cLambda1HbO*log(hdChannel1)-settings.hd.cLambda2HbO*log(hdChannel2));
            Hb=(tmpA0Hb+settings.hd.cLambda1Hb*log(hdChannel1)-settings.hd.cLambda2Hb*log(hdChannel2));
            clearvars tmp* i*
            toc;

            HbO(isnan(HbO)) = 0;
            Hb(isnan(Hb)) = 0;
            HbO = f_smooth2d_WF(HbO,2);
            Hb = f_smooth2d_WF(Hb,2);
            
            if commonSettings.saveHDF5
                fprintf('\nSaving Hb and HbO data (h5)...')
                h5create([commonFolders.processed,filesep,'run',num2str(dataIn(folderCounter).runnum,'%04.0f'),'.h5'],'/hemodynamics/Hb',size(Hb),'Deflate',commonSettings.compression,'Chunksize',commonSettings.chunkSize,'Datatype','double')
                h5create([commonFolders.processed,filesep,'run',num2str(dataIn(folderCounter).runnum,'%04.0f'),'.h5'],'/hemodynamics/HbO',size(HbO),'Deflate',commonSettings.compression,'Chunksize',commonSettings.chunkSize,'Datatype','double')
                h5write([commonFolders.processed,filesep,'run',num2str(dataIn(folderCounter).runnum,'%04.0f'),'.h5'],'/hemodynamics/Hb',Hb)
                h5write([commonFolders.processed,filesep,'run',num2str(dataIn(folderCounter).runnum,'%04.0f'),'.h5'],'/hemodynamics/HbO',HbO)
                fprintf('done.\n');toc
            end
        end
        
        %% process gfp

        if ~isempty(settings.gfpChannel)
            fprintf('\nProcess green fluorescence...')

            tmpWL = [470 515];
            tmpExtinction = f_GetExtinctions_WF(tmpWL);     % in cm
            tmpPathEx = f_pathlengths_WF(tmpWL(1),0.4)/2;       % pathlengths returns in cm
            tmpPathEm = f_pathlengths_WF(tmpWL(2),0.4)/2;
            tmpMuaEx = (tmpExtinction(1,1).*HbO) + (tmpExtinction(1,2).*Hb);
            tmpMuaEm = (tmpExtinction(2,1).*HbO) + (tmpExtinction(2,2).*Hb);
            gfp_norm = gfpChannel./mean(gfpChannel,3);
            gfp_norm_HD = (gfp_norm).*exp((tmpMuaEx.*tmpPathEx + tmpMuaEm.*tmpPathEm))-1;
            gfp_norm = gfp_norm - 1;

            gfp_norm(isnan(gfp_norm)) = 0;
            gfp_norm_HD(isnan(gfp_norm_HD)) = 0;
            gfp_norm = f_smooth2d_WF(gfp_norm,2);
            gfp_norm_HD = f_smooth2d_WF(gfp_norm_HD,2);

            clearvars tmp* i*
            toc;
            if commonSettings.saveHDF5
                fprintf('\nSaving green fluorescence data (h5)...')
                h5create([commonFolders.processed,filesep,'run',num2str(dataIn(folderCounter).runnum,'%04.0f'),'.h5'],'/gfp/norm',size(gfp_norm),'Deflate',commonSettings.compression,'Chunksize',commonSettings.chunkSize,'Datatype','double');
                h5write([commonFolders.processed,filesep,'run',num2str(dataIn(folderCounter).runnum,'%04.0f'),'.h5'],'/gfp/norm',gfp_norm)
                h5create([commonFolders.processed,filesep,'run',num2str(dataIn(folderCounter).runnum,'%04.0f'),'.h5'],'/gfp/normHD',size(gfp_norm_HD),'Deflate',commonSettings.compression,'Chunksize',commonSettings.chunkSize,'Datatype','double');
                h5write([commonFolders.processed,filesep,'run',num2str(dataIn(folderCounter).runnum,'%04.0f'),'.h5'],'/gfp/normHD',gfp_norm_HD)
                fprintf('done.\n');toc
            end
        end
        toc;
        
        %% process rfp
        if ~isempty(settings.rfpChannel)
            %%
            fprintf('\nProcess red fluorescence...')
            %optional: add filter and other pp
            
            if ~isempty(settings.hdChannel1)
                tmpHbRed = hdChannel2./mean(hdChannel2,3);
                tmpHbGreen = hdChannel1./mean(hdChannel1,3);
            end
            rfp_norm = rfpChannel./mean(rfpChannel,3);
            if ~isempty(settings.hdChannel1)
                rfp_norm_HD = rfp_norm./(tmpHbRed.^0.8.*tmpHbGreen.^0.4) - 1;
            end
            rfp_norm = rfp_norm-1;

            rfp_norm(isnan(rfp_norm)) = 0;
            rfp_norm_HD(isnan(rfp_norm_HD)) = 0;
            rfp_norm = f_smooth2d_WF(rfp_norm,2);
            rfp_norm_HD = f_smooth2d_WF(rfp_norm_HD,2);

            clearvars tmp* i
            toc;
            if commonSettings.saveHDF5
                fprintf('\nSaving red fluorescence data (h5)...')
                h5create([commonFolders.processed,filesep,'run',num2str(dataIn(folderCounter).runnum,'%04.0f'),'.h5'],'/rfp/norm',size(rfp_norm),'Deflate',commonSettings.compression,'Chunksize',commonSettings.chunkSize,'Datatype','double')
                h5write([commonFolders.processed,filesep,'run',num2str(dataIn(folderCounter).runnum,'%04.0f'),'.h5'],'/rfp/norm',rfp_norm)
                if ~isempty(settings.hdChannel1)
                    h5create([commonFolders.processed,filesep,'run',num2str(dataIn(folderCounter).runnum,'%04.0f'),'.h5'],'/rfp/normHD',size(rfp_norm_HD),'Deflate',commonSettings.compression,'Chunksize',commonSettings.chunkSize,'Datatype','double')
                    h5write([commonFolders.processed,filesep,'run',num2str(dataIn(folderCounter).runnum,'%04.0f'),'.h5'],'/rfp/normHD',rfp_norm_HD)
                end
                fprintf('done.\n');toc
            end
        end
    end
end

end