%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            f_load_ephys
% author - Brad Rauscher (created 2026)
% 
% Loads electrophysiology data acquired with the Intan Electrophyisology
% software.
% 
% INPUTS: f_load_ephys(mouse, date, run, _)
%   mouse: mouse name
%   date: date of acquisition
%   run: run number (int)
% 
% OPTIONAL INPUTS:
%   trim: trim outputs with digital input (default = 0)
%   skip: number of frames to skip to correct 'trim' (default = 1)
% 
% OUTPUTS:
%   ECoG: raw amplifier data (T x N)
%   t: time axis
%   imp: impedance values for each channel in MOhms (N)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ECoG, t, imp, fs] = f_load_ephys(mouse, date, run, varargin)

% parse inputs
p = inputParser;
addParameter(p, 'trim', false);
addParameter(p, 'skip', 1);

parse(p, varargin{:});

% organize folders
files = struct;
files.Root = '/projectnb/devorlab/bcraus/HRF/1P';
files.run = fullfile(files.Root, date, mouse);
files.rhd = fullfile(files.run, 'ephys');

for i = 1 : 2
    runs = dir(files.rhd);
    runs(1 : 2) = [];
    idx = contains({runs.name}, sprintf('Run%02i', run), IgnoreCase = true);
    files.rhd = fullfile(files.rhd,runs(idx).name);
end

rhd = f_read_rhd_2022(files.rhd);

% organize data
imp = [rhd.amplifier_channels.electrode_impedance_magnitude]' / 1e6;

trigger = find(diff(rhd.board_dig_in_data) == 1);

if p.Results.trim
    period = round(mean(diff(trigger(1 : p.Results.skip : end))));
else
    period = 0;
end

trigger = trigger(1) : trigger(end) + period;

ECoG = rhd.amplifier_data(:, trigger)';

fs = rhd.frequency_parameters.amplifier_sample_rate;

t = (1 : size(ECoG, 1))' / fs;

end