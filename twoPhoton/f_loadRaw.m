%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              f_loadRaw
% author - Brad Rauscher (created 2026)
% 
% Loads raw, motion-corrected data (data.bin) from suite2P output. Requires
% "ops.npy" file.
% 
% INPUTS: f_loadRaw(path)
%   path: path to parent directory of "data.bin" containing "ops.npy"
% 
% OUTPUTS:
%   data_raw: raw, motion-corrected data
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [chan1, chan2] = f_loadRaw(path)

contents = dir(path);
flag_1 = 0;
flag_2 = 0;

if ismember('data.bin',{contents.name})
    flag_1 = 1;
    fprintf('Found data.bin\n');
end

if ismember('data_chan2.bin',{contents.name})
    flag_2 = 1;
    fprintf('Found data_chan2.bin\n');
end

if ismember('ops.npy',{contents.name})
    fprintf('Found ops.npy');
end

folder = fileparts(mfilename('fullpath'));

pyenv('Version','/projectnb/devorlab/bcraus/envs/.venv/bin/python');
np = py.importlib.import_module('numpy');
ops = np.load(fullfile(path,'ops.npy'), allow_pickle=true);

ops = ops.item();

Ly = double(ops{'Ly'});
Lx = double(ops{'Lx'});

% load chan 1
if flag_1
    fid = fopen(fullfile(path,'data.bin'), 'r');
    data_raw = fread(fid, 'int16');
    fclose(fid);
    
    data_raw = reshape(data_raw, Lx, Ly, []);
    chan1 = permute(data_raw, [2, 1, 3]);
end

% load chan 2
if flag_2
    fid = fopen(fullfile(path,'data_chan2.bin'), 'r');
    data_raw = fread(fid, 'int16');
    fclose(fid);
    
    data_raw = reshape(data_raw, Lx, Ly, []);
    chan2 = permute(data_raw, [2, 1, 3]);
end

end