function [rawData] = f_loadRaw1P(Mouse,Date,Run,project,user,nChannels,frames)

Root_folder = '/projectnb/devorlab';
Root_folder = fullfile(Root_folder,user,project,Date,Mouse);

tmpSettings.mp = 1;
tmpSettings.cores= 4;
tmpSettings.doNorm = 0;
tmpSettings.doCat=0;
tmpSettings.nImport = 0;
tmpSettings.nChannels = nChannels;
if nargin > 6
    tmpSettings.frames = frames;
end
rawData = f_AndorDATImporter(fullfile(Root_folder,'onephoton',sprintf('Run%02i',Run)),tmpSettings);
