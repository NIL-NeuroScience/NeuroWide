function behCam = f_loadBeh(path)

files = dir([path,'/*.tiff']);
files = sort(fullfile({files.folder},{files.name}))';

frame_num = regexp(files,'\d+','match');
frame_num = cellfun(@(x) str2double(x{end}), frame_num);

[~,order] = sort(frame_num);

files = files(order);

N = numel(files);

template = Tiff(files{1},'r');
template = im2uint8(read(template));

behCam = zeros([size(template),N]);

for i = 1:N
    t = Tiff(files{i},'r');
    behCam(:,:,i) = im2uint8(read(t));
end

end