function data = f_load_compressed_mp4(path)

folder = f_path();

pyenv('Version',fullfile(folder,'.venv/bin/python'));
imageio = py.importlib.import_module('imageio.v2');
np = py.importlib.import_module('numpy');

reader = imageio.get_reader(path, 'ffmpeg');

meta = reader.get_meta_data();
T = double(reader.count_frames());

dim = double(meta{'size'});
H = dim(2);
W = dim(1);

data = zeros(H, W, T);

for i = 0:T-1
    frame = double(reader.get_data(py.int(i)));
    data(:,:,i+1) = frame(:,:,1);
end
reader.close()

end