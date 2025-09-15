function [nwb,data,Signals,brain_mask,allen,acq_settings,triggers,behCam] = f_loadNWB(date,mouse,run,user,root)

%% gather and organize inputs

pSet = struct;
pSet.date = date;
pSet.mouse = mouse;
pSet.run = run;

if nargin < 3; error('!!! Not enough inputs !!!'); end
if nargin < 4; pSet.user = 'bcraus'; else; pSet.user = user; end
if nargin < 5; pSet.root = 'HRF/1P'; else; pSet.root = root; end

%% setup files

files = struct;
files.root = fullfile('/projectnb/devorlab',pSet.user,pSet.root);
files.nwb = fullfile(files.root,pSet.date,pSet.mouse,'nwb',sprintf('Run%02i.nwb',pSet.run));
files.nwbPackage = '/projectnb/devorlab/bcraus/AnalysisCode/new_processing/matnwb';

if ~f_checkFile_WF(files.nwb)
    error('!!! NWB file not found !!!');
end

addpath(files.nwbPackage);

nwb = nwbRead(files.nwb,'ignorecache');

%% load videos

data = struct; 

data.rfp = nwb.acquisition.get('rfp').data.load();
data.rfp_HD = nwb.acquisition.get('rfp_HD').data.load();
data.gfp = nwb.acquisition.get('gfp').data.load();
data.gfp_HD = nwb.acquisition.get('gfp_HD').data.load();
data.HbO = nwb.acquisition.get('HbO').data.load();
data.HbR = nwb.acquisition.get('HbR').data.load();
data.HbT = data.HbO+data.HbR;

%% triggers

triggers = nwb.acquisition.get('DAQ_triggers').data.load();
triggers = array2table(triggers','VariableNames',split(nwb.acquisition.get('DAQ_triggers').description,';'));

%%
cIdx = struct;
cIdx.c470 = find(triggers.LED470Signal);
cIdx.c565 = find(triggers.LED565Signal);
cIdx.c525 = find(triggers.LED525Signal);
cIdx.c625 = find(triggers.LED625Signal);

channels = {'470','565','525','625'};
channel_order = zeros(1,4);
for i = 1:4
    if ~isempty(cIdx.(['c' channels{i}]))
        channel_order(i) = cIdx.(['c' channels{i}])(1);
    else
        channel_order(i) = 0;
    end
end

channels = channels(find(channel_order));
channel_order(channel_order==0) = [];
[~,channel_order] = sort(channel_order);
channels = channels(channel_order);

%%

acq_settings = struct;
acq_settings.channels = keys(nwb.general_optophysiology);
first_channel = acq_settings.channels{1};

acq_settings.channels = cellfun(@(s) s(1:3),acq_settings.channels,'UniformOutput',false);
if ~isequal(sort(channels),sort(acq_settings.channels))
    error('!!! Channels do not match (triggers vs. nwb) !!!');
end

acq_settings.fs = nwb.general_optophysiology.get(first_channel).imaging_rate;
acq_settings.DAQ_frequency = nwb.acquisition.get('DAQ_triggers').starting_time_rate;

acq_channels = keys(nwb.acquisition);
acq_channels(strcmp(acq_channels,'BehaviorVideo')) = [];
acq_channels(strcmp(acq_channels,'DAQ_triggers')) = [];
acq_channels(contains(acq_channels,'HD')) = [];
for i = 1:numel(acq_channels)
    acq_channels_planes{i} = nwb.acquisition.get(acq_channels{i}).imaging_plane.path;
    acq_channels_power(i) = nwb.acquisition.get(acq_channels{i}).power;
    acq_channels_exposure(i) = nwb.acquisition.get(acq_channels{i}).exposure_time;
end
acq_channels_planes = cellfun(@(s) s(end-11:end-9),acq_channels_planes,'UniformOutput',false);
[~,idx] = ismember(channels,acq_channels_planes);
acq_settings.LED_Power = acq_channels_power(idx);
acq_settings.LED_Exposure = acq_channels_exposure(idx);

%% behavior

beh_keys = keys(nwb.processing.get('behavior').nwbdatainterface.get('BehavioralTimeSeries').timeseries);
Signals = struct;
if find(ismember(beh_keys,'pupil_diameter'))
    Signals.pupil = nwb.processing.get('behavior').nwbdatainterface.get('BehavioralTimeSeries').timeseries.get('pupil_diameter').data.load();
end
if find(ismember(beh_keys,'whisker_long'))
    Signals.whisker_long = nwb.processing.get('behavior').nwbdatainterface.get('BehavioralTimeSeries').timeseries.get('whisker_long').data.load();
end
if find(ismember(beh_keys,'whisker_pad'))
    Signals.whisker_pad = nwb.processing.get('behavior').nwbdatainterface.get('BehavioralTimeSeries').timeseries.get('whisker_pad').data.load();
end
if find(ismember(beh_keys,'accelerometer'))
    Signals.accelerometer = nwb.processing.get('behavior').nwbdatainterface.get('BehavioralTimeSeries').timeseries.get('accelerometer').data.load();
end

%% allen ROIs

allen_keys = keys(nwb.processing.get('ophys').nwbdatainterface);
allen_keys(strcmp(allen_keys,'ImageSegmentation')) = [];

for i = 1:numel(allen_keys)
    Signals.allen.(allen_keys{i}) = nwb.processing.get('ophys').nwbdatainterface.get(allen_keys{i}).data.load()';
end

%% extract masks

dim = nwb.acquisition.get(acq_channels{1}).dimension.load();

allen_nwb = nwb.processing.get('ophys').nwbdatainterface.get('ImageSegmentation').planesegmentation.get('Allen_segmentation').pixel_mask.data.load();
allen_idx = nwb.processing.get('ophys').nwbdatainterface.get('ImageSegmentation').planesegmentation.get('Allen_segmentation').pixel_mask_index.data.load();
allen_idx = [0;allen_idx];

allen_labels = split(nwb.processing.get('ophys').nwbdatainterface.get('ImageSegmentation').planesegmentation.get('Allen_segmentation').pixel_mask.description,';');

mask = NaN([dim' 12]);

for mask_idx = 1:12
    x = allen_nwb.x(allen_idx(mask_idx)+1:allen_idx(mask_idx+1));
    y = allen_nwb.y(allen_idx(mask_idx)+1:allen_idx(mask_idx+1));
    for i = 1:numel(x)
        mask(y(i),x(i),mask_idx) = 1;
    end
end

allen = struct;
allen.masks = mask;
allen.labels = allen_labels;

brain_mask = NaN(dim');
BM = nwb.processing.get('ophys').nwbdatainterface.get('ImageSegmentation').planesegmentation.get('brain_mask').pixel_mask.data.load();
x = BM.x;
y = BM.y;
for i = 1:numel(x)
    brain_mask(y(i),x(i)) = 1;
end

%% extract behavior video

if nargout > 7
    behCam = nwb.acquisition.get('BehaviorVideo').data.load();
end

end