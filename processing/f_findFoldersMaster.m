function folderList=f_findFoldersMaster(folderIn,type)
%Extracts a list of folders of a specific type from folder structure.
%Filelist needs to be a global variable to be accessible through all levels of recursion.
global fileList;fileList={};f_findFolderByType(folderIn,type)
if ~isempty(fileList)
    folderList=cell2table(fileList,'VariableNames',{'type','parentFolder','name','runnum'});
else
    folderList=[];
end
end