function [exist] = f_checkFile(file)

exist = false;

folders = split(file,'/');
parent = join(folders(1:end-1),'/');
file = folders(end);

tmpDir = dir(string(parent));
tmpDir = find(strcmp({tmpDir.name},string(file))==1,1);

if ~isempty(tmpDir)
    exist = ~exist;
end