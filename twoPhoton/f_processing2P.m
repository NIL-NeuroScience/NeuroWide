%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%                     Two photon data processing
% 
% Initial processing script for two photon imaging data. 
% 
% Authors: Brad Rauscher (April 2026)
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

function f_processing2P(Mouse, Date, varargin)
%%
Mouse = 'Rbp4_139';
Date = '26-04-13';
varargin = {};

%% parse inputs
p = inputParser;
addParameter(p,'runs',[]);
% addParameter(p,'smooth',2);
% addParameter(p,'compression',0);
% addParameter(p,'rotation',0);
addParameter(p,'load_dir','bcraus/HRF/2P');
addParameter(p,'save_dir','bcraus/HRF/2P');

parse(p,varargin{:});

%% organize files
path = fullfile('/projectnb/devorlab', ...
    p.Results.load_dir, ...
    Date, ...
    Mouse);
path_data = fullfile(path, 'twophoton');
path_behavior = fullfile(path, 'camera');
path_ephys = fullfile(path, 'ephys');

% find runs
N = numel(p.Results.runs);
runs = [];
if ~N
    run_names = dir(path_data);
    run_names = {run_names.name};

    for i = 1:numel(run_names)
        if contains(run_names{i}, 'Run')
            run_num_str = strsplit(run_names{i}, 'Run');
            runs = [runs, str2double(run_num_str{2}(1:2))];
        end
    end
else
    runs = p.Results.runs;
end

if exist(path_ephys)
    ephys_contents = dir(path_ephys);
end

dataIn = struct(runnum=[], twophoton=[], behavior=[], ephys=[], settings=[], template=[]);

runs = sort(runs);
for i = 1:numel(runs)
    dataIn(i).runnum = runs(i);
    
    % organize twophoton files
    dataIn(i).twophoton.runnum = runs(i);
    dataIn(i).twophoton.name = run_names{contains(run_names, sprintf('Run%02i', runs(i)))};
    dataIn(i).twophoton.folder = path_data;

    % organize behavior videos
    flag_beh = exist(fullfile(path_behavior, sprintf('Run%02i.mp4', runs(i))));
    if flag_beh
        dataIn(i).behavior.runnum = runs(i);
        dataIn(i).behavior.name = sprintf('Run%02i.mp4', runs(i));
        dataIn(i).behavior.folder = path_behavior;
    else
        dataIn(i).behavior = [];
    end

    % organize ephys files
    flag_beh = exist(fullfile(path_behavior, sprintf('Run%02i.mp4', runs(i))));
    if flag_beh
        dataIn(i).behavior.runnum = runs(i);
        dataIn(i).behavior.name = sprintf('Run%02i.mp4', runs(i));
        dataIn(i).behavior.folder = path_behavior;
    else
        dataIn(i).behavior = [];
    end
end


%%
contents = dir(path_data);
path_run = contents(contains({contents.name}, sprintf('Run%02i', Run)))

end