function [folders, files] = f_getFolderContent(in,filetype)
tmpList=dir(in);
for i=1:size(tmpList,1)
    if tmpList(i).isdir && tmpList(i).name(1) ~= '.'
        tmpFolders{i}=tmpList(i).name;
    elseif ~tmpList(i).isdir && tmpList(i).name(1) ~= '.'
        if strcmp(filetype,'*')
            tmpFiles{i}=tmpList(i).name;
        else
            if strcmp(tmpList(i).name(end-2:end),filetype)
                tmpFiles{i,1}=tmpList(i).name;
            end
        end
    end
end
if exist('tmpFolders','var')
    folders=tmpFolders(~cellfun('isempty',tmpFolders));
else
    folders={};
end
if exist('tmpFiles','var')
    files=tmpFiles(~cellfun('isempty',tmpFiles));
else
    files={};
end
