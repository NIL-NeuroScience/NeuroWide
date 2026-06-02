function rawData = loadRaw1P(path, nChannels, frames)

tmpSettings.mp = 1;
tmpSettings.cores= 4;
tmpSettings.doNorm = 0;
tmpSettings.doCat=0;
tmpSettings.nImport = 0;
tmpSettings.nChannels = nChannels;
if nargin > 6
    tmpSettings.frames = frames;
end

rawData = NeuroWide.io.andorDATImporter(path, tmpSettings);
