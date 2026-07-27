
function vid = loadTifStack(path)

% collect metadata
meta = imfinfo(path);
n_frames = numel(meta);
H = meta(1).Height;
W = meta(1).Width;
bit_depth = meta(1).BitDepth;

% set proper datatype
if bit_depth == 8
    dtype = 'uint8';
elseif bit_depth == 16
    dtype = 'uint16';
end

% initialize vid and tiff
vid = zeros(H, W, n_frames, dtype);
t = Tiff(path, 'r');

for i = 1:n_frames
    t.setDirectory(i);
    vid(:,:,i) = t.read();
end

t.close();

end