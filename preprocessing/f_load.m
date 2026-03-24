%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%                               f_load.m
% 
% Secondary processing script for widefield fluorescence and intrinsic
% imaging data. Creates brain and vessel masks, estimates pupil diameter,
% whisking, and movement, and performs preliminary analysis between the
% optical signals.
% 
% INPUTS: f_load(Mouse,Date,Run,idx,_)
%   Mouse: mouse ID
%   Date: session date
%   Runs: list of run to load
%   idx: index of which signals to load [rfp, rfp_HD, gfp, gfp_HD, HbO,
%       HbR] (loading HbO and HbR will automatically calculate HbT)
% 
% EXTRA PARAMETERS:
%   load_dir: directory to load from
%   behCam: nothing (to be added in the future)
%   range: frame indices to load
%   cores: number of cores to use to load data
% 
% OUTPUTS:
%   data: contains processed widefield videos
%   Signals: contains behavioral measurements
%   brain_mask: mask of the window exposure
%   parcellation: structure containing individual allen atlas regions
%   vessel_mask: mask of large pial vessels
%   dataIn: info about runs
%   settings: acquisition settings (contains sampling rate)
%   files: associated files
%   digitalInput: acquisition timing data
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [data,Signals,brain_mask,parcellation,vessel_mask,dataIn,...
    settings,files,digitalInput] = f_load(Mouse,Date,Run,idx,varargin)

%% parse inputs
p = inputParser;
addParameter(p,'load_dir','/projectnb/devorlab/bcraus/HRF/1P');
addParameter(p,'behCam',false);
addParameter(p,'range',[]);
addParameter(p,'cores',4);
parse(p,varargin{:});

%% organize files
files.Root_Folder = p.Results.load_dir;
files.Root_Folder = fullfile(files.Root_Folder,Date,Mouse);
files.h5 = fullfile(files.Root_Folder,'processed',sprintf('run%04i.h5',Run));
files.dataIn = fullfile(files.Root_Folder,'dataIn.mat');
files.Signals = fullfile(files.Root_Folder,'DataAnalysis',sprintf('Run%02i',Run),'Signals.mat');
files.ROIs = fullfile(files.Root_Folder,'DataAnalysis','ROIs.mat');
files.triggers = fullfile(files.Root_Folder,'Triggers',sprintf('Run%03i.mat',Run));

%% load masks
if f_checkFile(files.ROIs)
    ROIs = load(files.ROIs);
else
    warning("'ROIs.mat' could not be found!")
    ROIs = struct;
end

if isfield(ROIs,'brain_mask')
    brain_mask = ROIs.brain_mask;
else
    warning("'brain_mask' is not a fieldname in ROIs! Returning []")
    brain_mask = [];
end

if isfield(ROIs,'vessel_mask')
    vessel_mask = ROIs.vessel_mask;
else
    warning("'vessel_mask' is not a fieldname in ROIs! Returning []")
    vessel_mask = [];
end

if isfield(ROIs,'parcellation')
    parcellation = ROIs.parcellation;
else
    warning("'parcellation' is not a fieldname in ROIs! Returning []")
    parcellation = [];
end

%% load dataIn

if f_checkFile(files.dataIn)
    dataIn = load(files.dataIn);dataIn = dataIn.dataIn;
else
    warning("'dataIn.mat' could not be found! Returning []")
    dataIn = [];
end

%% load signals

if f_checkFile(files.Signals)
    Signals = load(files.Signals);Signals = Signals.Signals;
else
    warning("'Signals.mat' could not be found! Returning []")
    Signals = [];
end

%% load digitalInput and settings

if f_checkFile(files.triggers)
    triggers = load(files.triggers);
else
    warning("Trigger file could not be found!")
    triggers = struct;
end

if isfield(triggers,'settings')
    settings = triggers.settings;
else
    warning("'settings' is not a fieldname in triggers! Returning []")
    settings = [];
end

if isfield(triggers,'digitalInput')
    digitalInput = triggers.digitalInput;
else
    warning("'digitalInput' is not a fieldname in triggers! Returning []")
    digitalInput = [];
end

%%

fIdx = find(idx);

dataOrder = {'/rfp/norm';'/rfp/normHD';'/gfp/norm';'/gfp/normHD';'/hemodynamics/HbO';'/hemodynamics/Hb'};
dataName = {'rfp';'rfp_HD';'gfp';'gfp_HD';'HbO';'HbR';'HbT'};

dataOrder = {dataOrder{fIdx}};
dataName = {dataName{fIdx}};

N = numel(dataOrder);

data = struct;
tmp = cell(N,1);

h5 = files.h5;

if isempty(gcp('nocreate'))
    parpool(p.Results.cores);
end

if ~isempty(p.Results.range)
    info = h5info(files.h5,dataOrder{1});
    dims = info.Dataspace.Size;
    numFrames = numel(p.Results.range);

    start = [1,1,p.Results.range(1)];
    count = [dims(1),dims(2),numFrames];

    parfor i = 1:N
        tmp{i} = h5read(h5,dataOrder{i},start,count);
    end
else
    parfor i = 1:N
        tmp{i} = h5read(h5,dataOrder{i});
    end
end


%%

for i = 1:N
    if fIdx(i) < 5
        data.(dataName{i}) = 100*tmp{i};
    else
        data.(dataName{i}) = 1e6*tmp{i};
    end
end

if idx(5) && idx(6)
    data.HbT = data.HbO+data.HbR;
end

end