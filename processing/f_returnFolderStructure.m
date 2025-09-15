function [folders,files]=f_returnFolderStructure(folderIn,level)
%Replicates a folder/file structure to a certain depth (level) into a struct array.
%Contains classification of folder types by their name and extraxts run number
%23-05-30 added onephotonSona and onephotonZyla to the list for 1photon
%folder to make distinction afterwards
[tmpFolders,tmpFiles]=f_getFolderContent(folderIn,'*');

if ~isempty(tmpFolders)
    for iFolder=1:size(tmpFolders,2)
        folders(iFolder).parentFolder=folderIn;
        folders(iFolder).name=tmpFolders{iFolder};
        %Identifies runnum in folder name
        if contains(folders(iFolder).name,'run','IgnoreCase',true)
            try
                tmprun=textscan(folders(iFolder).name(regexp(folders(iFolder).name,'[Rr]un')+3:end),'%d%s');
                folders(iFolder).runnum=double(tmprun{1});
            catch
                folders(iFolder).runnum=nan;
            end
        else
            folders(iFolder).runnum=nan;
        end
        %Categorization of folder types
        if any(strcmpi(folders(iFolder).name,{'2photon','twophoton'}))
            folders(iFolder).type='2photon';
        elseif any(contains(folders(iFolder).name,{'ZSeries','TSeries','SingleImage','PointScan'}))
            folders(iFolder).type='2photonRaw';
        elseif any(contains(folders(iFolder).name,{'ephys','laminar','intan'}))
            folders(iFolder).type='ephys';
        elseif any(strcmpi(folders(iFolder).name,{'1photon','singlephoton','micromanager','solis','onephoton','onephotonSona','onephotonZyla'}))
            folders(iFolder).type='1photon';
        elseif any(strcmpi(folders(iFolder).name,{'trigger','triggers','daq'}))
            folders(iFolder).type='daq';
        elseif any(strcmpi(folders(iFolder).name,{'pylon','webcam'}))
            folders(iFolder).type='behaviorCam';
        elseif any(strcmpi(folders(iFolder).name,{'processed','processedEP','matlabdata'}))
            folders(iFolder).type='processed';
        elseif any(strcmpi(folders(iFolder).name,{'images','imagesEP','plots','videos'}))
            folders(iFolder).type='imageOutput';
        else
            folders(iFolder).type=[];
        end
        %adds recursion
        if level>1
            [folders(iFolder).folders,folders(iFolder).files]=f_returnFolderStructure([folders(iFolder).parentFolder,filesep,folders(iFolder).name],level-1);
        end
    end
else
    folders=[];
end

if ~isempty(tmpFiles)
    files=cell(size(tmpFiles,2),3);
    for iFile=1:size(tmpFiles,2)
        %filename
        files{iFile,1}=tmpFiles{iFile};
        %filetype
        tmpDots=strfind(tmpFiles{iFile},'.');
        if ~isempty(tmpDots)
            files{iFile,2}=tmpFiles{iFile}(tmpDots(end)+1:end);
        else
            files{iFile,2}='';
        end
        %runnum
        if contains(tmpFiles{iFile},'run','IgnoreCase',true)
            tmprun=textscan(tmpFiles{iFile}(regexp(tmpFiles{iFile},'[Rr]un')+3:end),'%d%s');
            files{iFile,3}=double(tmprun{1});
        else
            files{iFile,3}=nan;
        end
    end
else
    files=[];
end

end