function rawImageChannels=f_AndorDATImporter(folderIn,settings)

%ANDOR SOLIS SPOOL.DAT IMPORT INTO MAtLAB
% Martin Thunemann, 2023-01-23, Version 2
% Brad Rauscher, 2023-10-30, Version 3
%   - Added correction for metadata appended to end of every image
% 
% f_AndorDATImporter imports time series data from Andor Solis software
% into MATLAB.
% 
% f_AndorDATImporter() will process data in current directory using
% standard settings and shows a montage of the DF/F and DR/R time series.
% f_AndorDATImporter(folderIn) will process data within folderIn using standard settings.
% f_AndorDATImporter(folderIn,settings) will process data within folderIn using user-defined settings.
% 
% INPUT ARGUMENTS (optional):
% folderIn (char), data folder to be imported.
% settings (struct), containing import settings (see standard settings)
% 
% OUTPUT ARGUMENTS (0-4):
% rawImage, imported image data (W x H x T).
% metadata, struct derived from image import.
% rawImageChannels, imported image data sorted by channels (W x H x T X C)
% combinedImageNorm, imported image data normalized by temporal average per
%   channel and combined into single image.
% 
% The input (i.e., data) folder needs to contain 
% -     XXXXXXXXXXspool.dat files
% -     acquisitionmetadata.ini
% -     Spooled files.sifx
% 

%% Standard settings 
stdSettings.mp = 1;        % run with parallel processing
stdSettings.cores= 4;      % number of CPU cores
stdSettings.doNorm = 1;    % optional: run normalization
stdSettings.nChannels = 4; % optional: number of channels per cycle
stdSettings.nImport=0;      %number of images to import (0=all, per channel)
stdSettings.doCat=1;       % optional: create image montage
% Define execution parameters

if ~ischar(folderIn) && ~isstruct(settings)
    fprintf('Error! Please check input variable(s)!'); error('Error!\nPlease check input variable(s)!');
end

files = struct;
files.meta = fullfile(folderIn,'acquisitionmetadata.ini');
files.spooled_files = fullfile(folderIn,'Spooled files.sifx');

%% Get started
tic;
fprintf('\n\nReading %s... \n',folderIn);
if exist(files.meta,'file')~=2 || exist(files.spooled_files,'file')~=2
    error('Error! acquisitionmetadata.ini or Spooled files.sifx missing!');
end

%% Starting Multiprocessor environment
if settings.mp
    poolobj = gcp('nocreate');
    if isempty(poolobj)
        parpool('local',settings.cores);
    end
end

% Read settings from acquisitionmetadata.ini
metadata=ini2struct(files.meta);

% Add current folder into metadata
metadata.folder=cd;

% Read some metadata from "spooled files.sifx"
[metadata.totalImages_sifx,metadata.acqusitionSettings] = f_sifx2struct(files.spooled_files);

% metadata.aoistride is given in Bytes; the metadata.aoistrideFactor depends on metadata.pixelencoding
% found on: https://andor3.readthedocs.io/en/latest/_modules/andor3/utils.html#decode_image_data

switch metadata.pixelencoding 
    case 'Mono32'
        metadata.aoistrideFactor=4;
        metadata.conversionSetting='uint32=>uint32';
        metadata.rawPixelFormat='uint32';
    case 'Mono12'
        metadata.aoistrideFactor=2;
        metadata.conversionSetting='uint16=>uint16';
        metadata.rawPixelFormat='uint16';
    case 'Mono16'
        metadata.aoistrideFactor=2;
        metadata.conversionSetting='uint16=>uint16';
        metadata.rawPixelFormat='uint16';
    otherwise
        fprintf('\nError: unkown pixel encoding\n');
        error('Error! Unkown pixel encoding!');
end
%metadata.aoiDataSize= metadata.aoistride/metadata.aoistrideFactor * metadata.imagesizebytes/metadata.aoistride * metadata.imagesperfile;
metadata.aoiDataSize= 1/metadata.aoistrideFactor * metadata.imagesizebytes * metadata.imagesperfile;

%% Find all .dat files and put file list into correct order
tmpFileList=dir(fullfile(folderIn,'*spool.dat'));
tmpFileList={tmpFileList.name};
fileImport.totalFiles=size(tmpFileList,2);
fileImport.totalImages=size(tmpFileList,2)*metadata.imagesperfile;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Taken from Columbia SCAPE_DataLoader_v1.m line 48ff %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmpFileStr=zeros(fileImport.totalFiles-1,10);
iSpool=0;
for iFile = 1:fileImport.totalFiles-1
    iSpool=iSpool+1;
    tmpFileNumber = iFile;
    for iExp = 1:10
        tmpFileStr(iFile, iExp) = mod(tmpFileNumber, 10^iExp)/(10^(iExp-1));
        tmpFileNumber = tmpFileNumber-mod(tmpFileNumber, 10^iExp);
    end
    tmpName = mat2str(tmpFileStr(iFile, :));
    tmpName = tmpName(2:end-1);
    tmpName = tmpName(tmpName ~= ' ');
    tmpName = fullfile(folderIn,[tmpName 'spool.dat']);
    fileImport.fileNames{iSpool} = tmpName;
end
fileImport.fileNames = [{fullfile(folderIn,'0000000000spool.dat')} fileImport.fileNames];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\nImages acquired: %d; Images stored: %d.',metadata.totalImages_sifx,fileImport.totalImages);
fprintf('\nImage dimensions: W %d x H %d x T %d (%s depth at %s encoding).', metadata.aoiwidth,metadata.aoiheight,metadata.totalImages_sifx,metadata.rawPixelFormat,metadata.pixelencoding)
fprintf('\nImage size: ~%0.2f MB.', metadata.imagesizebytes*fileImport.totalImages/(1024^2))
clearvars i* tmp*
%% Estimating how many images to be imported (new 1/26/23)
if settings.nImport==0
    fileImport.imagesRequested=metadata.totalImages_sifx;
    fileImport.filesRequested=fileImport.totalFiles;
    fileImport.imagesImported=fileImport.totalImages;
else
    fileImport.imagesRequested=settings.nImport*settings.nChannels;
    fileImport.filesRequested=ceil(fileImport.imagesRequested/metadata.imagesperfile);
    fileImport.imagesImported=fileImport.filesRequested*metadata.imagesperfile;
    if fileImport.filesRequested>=fileImport.totalFiles
        fileImport.imagesRequested=metadata.totalImages_sifx;
        fileImport.filesRequested=fileImport.totalFiles;
        fileImport.imagesImported=fileImport.totalImages;
        fprintf('\nImporting entire image series (%i images)',fileImport.imagesRequested)
    else
        fprintf('\nImporting partial image series (%i images)',fileImport.imagesRequested)
    end
end

%% Performing import

removeRows = metadata.aoiheight*metadata.aoiwidth*2;
removeRows = metadata.imagesizebytes-removeRows;

fprintf('\n\nImporting data ');
rawImage=zeros(fileImport.filesRequested,metadata.aoiwidth,metadata.aoiheight,metadata.imagesperfile,metadata.rawPixelFormat);
tmpFileHandler=zeros(fileImport.filesRequested,1);
if settings.mp %#ok<*PFBNS>
    fprintf('using %d CPUs...',settings.cores);
    parfor (iFile=1:fileImport.filesRequested,4)
        if removeRows == 40 || removeRows == 42
            tmpFileHandler(iFile)=fopen(fileImport.fileNames{iFile},'r');
            tmpIn=fread(tmpFileHandler(iFile),metadata.imagesizebytes*metadata.imagesperfile,metadata.conversionSetting,'l');
            metaInd = 0:(metadata.aoiDataSize/metadata.imagesperfile):metadata.aoiDataSize;metaInd(1) = [];
            tmpMetaInd = ones(metadata.aoiDataSize,1);for i = metaInd;tmpMetaInd((i-removeRows/2):(i-1)) = 0;end % use 20 for uint16 and 10 for uint32
            tmpRaw=reshape(tmpIn(logical(tmpMetaInd)),[metadata.aoistride/metadata.aoistrideFactor,(metadata.imagesizebytes-removeRows)/metadata.aoistride,metadata.imagesperfile]);
            rawImage(iFile,:,:,:)=tmpRaw(1:metadata.aoiwidth,1:metadata.aoiheight,:);
            fclose(tmpFileHandler(iFile));
        else
            tmpFileHandler(iFile)=fopen(fileImport.fileNames{iFile},'r');
            tmpIn=fread(tmpFileHandler(iFile),metadata.imagesizebytes*metadata.imagesperfile,metadata.conversionSetting,'l');
            tmpRaw=reshape(tmpIn(1:metadata.aoiDataSize),[metadata.aoistride/metadata.aoistrideFactor,metadata.imagesizebytes/metadata.aoistride,metadata.imagesperfile]);
            rawImage(iFile,:,:,:)=tmpRaw(1:metadata.aoiwidth,1:metadata.aoiheight,:);
            fclose(tmpFileHandler(iFile));
        end
    end
else
    fprintf('using single CPU...');
    for iFile=1:fileImport.filesRequested
        if removeRows == 40 || removeRows == 42
            tmpFileHandler(iFile)=fopen(fileImport.fileNames{iFile},'r');
            tmpIn=fread(tmpFileHandler(iFile),metadata.imagesizebytes*metadata.imagesperfile,metadata.conversionSetting,'l');
            metaInd = 0:(metadata.aoiDataSize/metadata.imagesperfile):metadata.aoiDataSize;metaInd(1) = [];
            tmpMetaInd = ones(metadata.aoiDataSize,1);for i = metaInd;tmpMetaInd((i-removeRows/2):(i-1)) = 0;end % use 20 for uint16 and 10 for uint32
            tmpRaw=reshape(tmpIn(logical(tmpMetaInd)),[metadata.aoistride/metadata.aoistrideFactor,(metadata.imagesizebytes-removeRows)/metadata.aoistride,metadata.imagesperfile]);
            rawImage(iFile,:,:,:)=tmpRaw(1:metadata.aoiwidth,1:metadata.aoiheight,:);
            fclose(tmpFileHandler(iFile));
        else
            tmpFileHandler(iFile)=fopen(fileImport.fileNames{iFile},'r');
            tmpIn=fread(tmpFileHandler(iFile),metadata.imagesizebytes*metadata.imagesperfile,metadata.conversionSetting,'l');
            tmpRaw=reshape(tmpIn(1:metadata.aoiDataSize),[metadata.aoistride/metadata.aoistrideFactor,metadata.imagesizebytes/metadata.aoistride,metadata.imagesperfile]);
            rawImage(iFile,:,:,:)=tmpRaw(1:metadata.aoiwidth,1:metadata.aoiheight,:);
            fclose(tmpFileHandler(iFile));
        end
    end
end
rawImage=permute(rawImage,[3,2,4,1]);
rawImage=reshape(rawImage,metadata.aoiheight,metadata.aoiwidth,fileImport.imagesImported);
% Remove 'empty' images at end of the series (mismatch dat/acquired images)
rawImage=rawImage(:,:,1:fileImport.imagesRequested);

tmpTime1 = toc;
fprintf('done after %0.2f s.\n',tmpTime1);
clearvars tmp* i*

%% Create normalized and/or concatenated image

fprintf('\nDividing image into %d channels...',settings.nChannels);
if mod(metadata.totalImages_sifx,settings.nChannels)~=0
    fprintf('Error! %d images cannot be divided in %d channels!\n',metadata.totalImages_sifx,settings.nChannels);
    frames_to_add = settings.nChannels - mod(metadata.totalImages_sifx,settings.nChannels);
    padding = rawImage(:,:,settings.nChannels+1-frames_to_add:settings.nChannels);
    rawImage = cat(3,padding,rawImage);
end
% Generate image with channel in dimension 4 (W x H x T x C)
rawImageChannels=zeros(metadata.aoiheight,metadata.aoiwidth,size(rawImage,3)/settings.nChannels,settings.nChannels,metadata.rawPixelFormat);
for iChannel=1:settings.nChannels
    rawImageChannels(:,:,:,iChannel)=rawImage(:,:,iChannel:settings.nChannels:end);
end
tmpTime2 = toc;
fprintf('done after %0.2f s.\n',tmpTime2);

rawImageChannels = double(rawImageChannels);

if settings.doNorm     
    % Normalize image across 3rd (T) dimension
    fprintf('\nPerforming normalization...');
    rawImageChannelsAvg=mean(rawImageChannels,3);
    rawImageChannelsNorm=100*(double(rawImageChannels)-rawImageChannelsAvg)./rawImageChannelsAvg;
    tmpTime3 = toc;
    fprintf('done after %0.2f s.\n',tmpTime3);

    % Create concatenated image
     fprintf('\nCreating image montage...');
    if settings.doCat
        if settings.nChannels==2
            combinedImageMean=cat(1,rawImageChannelsAvg(:,:,1),rawImageChannelsAvg(:,:,2));
            combinedImageNorm=cat(1,rawImageChannelsNorm(:,:,:,1),rawImageChannelsNorm(:,:,:,2));
        elseif settings.nChannels==3
            combinedImageMean=cat(1,rawImageChannelsAvg(:,:,1),rawImageChannelsAvg(:,:,2),rawImageChannelsAvg(:,:,3));
            combinedImageNorm=cat(1,rawImageChannelsNorm(:,:,:,1),rawImageChannelsNorm(:,:,:,2),rawImageChannelsNorm(:,:,:,3));
        elseif settings.nChannels==4
            combinedImageMean=cat(2,cat(1,rawImageChannelsAvg(:,:,1),rawImageChannelsAvg(:,:,3)),cat(1,rawImageChannelsAvg(:,:,2),rawImageChannelsAvg(:,:,4)));
            combinedImageNorm=cat(2,cat(1,rawImageChannelsNorm(:,:,:,1),rawImageChannelsNorm(:,:,:,3)),cat(1,rawImageChannelsNorm(:,:,:,2),rawImageChannelsNorm(:,:,:,4)));
        else
            combinedImageMean=[];
            combinedImageNorm=[];
        end
    tmpTime4 = toc;
    fprintf('done after %0.2f s.\n',tmpTime4);
    end
else
    combinedImageNorm=[];
end
clear tmp* i*
if nargout==0;figure();imagesc(combinedImageMean);axis image off;colorbar;colormap gray;implay(combinedImageNorm);return;end
end
%%
function out = ini2struct(FileName)
%modified from https://www.mathworks.com/matlabcentral/fileexchange/17177-ini2struct
out = [];
f = fopen(FileName,'r');
while ~feof(f)
    s = strtrim(fgetl(f));
    if isempty(s)
        continue;
    end
    if s(1)==';'
        continue;
    end
    if s(1)=='#'
        continue;
    end
    if any(contains(s,["[","]"]))
        continue;
    else
        [par,val] = strtok(s, '=');
        par = f_cleanvalue(par);
        val = f_cleanvalue(val);
        out.(lower(matlab.lang.makeValidName(par))) = val;
    end
end
fclose(f);
end

function res = f_cleanvalue(s)
%modified from https://www.mathworks.com/matlabcentral/fileexchange/17177-ini2struct
res = strtrim(s);
if size(res,1)>0 && strcmpi(res(1),'=')
    res(1)=[];
end
res = strtrim(res);
if ~isnan(str2double(res));res=str2double(res);end
end

function [nImages,acqusitionSettings] = f_sifx2struct(FileName)
nImages = [];
acqusitionSettings = [];
f = fopen(FileName,'r');
while ~feof(f)
    s = strtrim(fgetl(f));
    if isempty(s)
        continue;
    end
    if contains(s,'Pixel number')
        %disp(s)
        tmpOut=strsplit(s,' ');
        if size(tmpOut,2)==10
            nImages=str2double(tmpOut{7});
        end
    elseif contains(s,'<PreAmpGainText>')
        tmpOut=strsplit(s,{'<','>'});
        if size(tmpOut,2)==5
            acqusitionSettings.preAmpSetting=tmpOut{3};
        end
    elseif contains(s,'<ExtendedDynamicRange>')
        tmpOut=strsplit(s,{'<','>'});
        if size(tmpOut,2)==5
            acqusitionSettings.extendedDynamicRange=str2double(tmpOut{3});
        end
    end
end
fclose(f);
end