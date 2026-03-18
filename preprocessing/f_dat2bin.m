%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              f_dat2bin
% author - Brad Rauscher (created 2026)
% 
% Converts Andor Solid acquisition (.dat) files into a single binary 
% (.bin) file.
% 
% INPUTS: f_dat2bin(path, channels)
%   path: path to acuisition save folder
%   channels: channel order list
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f_dat2bin(path, channels)
%% organize files

meta_path = fullfile(path, 'meta.json');
bin_path = fullfile(path, 'data.bin');

%%
settings.mp = 1;
settings.cores= 4;
settings.doNorm = 0;
settings.doCat=0;
settings.nImport = 0;
settings.nChannels = numel(channels);
[rawData, dat_meta] = f_AndorDATImporter(path, settings);

% resave .bin
[h,w,T,c] = size(rawData);

%% save metadata .json

metadata = struct;
metadata.shape = [h, w, c, T];
metadata.dtype = dat_meta.rawPixelFormat;
metadata.order = 'F';
metadata.axes = {'H', 'W', 'C', 'T'};
metadata.units = 'intensity';
metadata.channel_order = channels;

jsonStr = jsonencode(metadata', PrettyPrint=true);
fid = fopen(meta_path, 'w');
fprintf(fid, '%s', jsonStr);
fclose(fid);

%% organize binary
binary = permute(rawData, [4,3,1,2]);
binary = reshape(binary, T*c, h, w);
binary = permute(binary,[2,3,1]);

binary = uint16(binary(:));

%% save binary

fid = fopen(bin_path, 'w');
fwrite(fid, binary, dat_meta.rawPixelFormat);
fclose(fid);

%% delete .dat files 
contents = dir(path);
contents(1:2) = [];

keep_idx = strcmp({contents.name}, 'data.bin') + strcmp({contents.name}, 'meta.json');
contents(logical(keep_idx)) = [];

for i = 1:numel(contents)
    if contents(i).isdir
        rmdir(fullfile(contents(i).folder, contents(i).name), 's')
    else
        delete(fullfile(contents(i).folder, contents(i).name));
    end
end

end