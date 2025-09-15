function f_findFolderByType(folderIn,type)
%Extracts a list of folders of a specific type from folder structure.
global fileList
for iFolder=1:size(folderIn,2)
    if ~isempty(folderIn(iFolder).type) && contains(folderIn(iFolder).type,type)
        fileList{end+1,1}=folderIn(iFolder).type;
        fileList{end,2}=folderIn(iFolder).parentFolder;
        fileList{end,3}=folderIn(iFolder).name;
        fileList{end,4}=folderIn(iFolder).runnum;
    end
    %adds recursion
    if isfield(folderIn(iFolder),'folders') && ~isempty(folderIn(iFolder).folders)
        f_findFolderByType(folderIn(iFolder).folders,type);
    end
end
end