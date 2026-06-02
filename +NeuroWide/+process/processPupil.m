%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
%    filenames: cell array of .tiff filepaths
%    behaviorROIs.pupil_mask: binary mask of eye
%    behaviorROIs.pupil_thrVal: threshold value to mask pupil vs. iris
%    behaviorROIs.eye_length: length (in pixels) horizontally across eye
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pupil] = processPupil(filenames, behaviorROIs)

tmpPupil = [];
tmpMask = behaviorROIs.pupil_mask;
tmpthrVal = behaviorROIs.pupil_thrVal;

if iscell(filenames)
    parfor i=1:numel(filenames)
        t = Tiff(filenames{i},'r');
        imageData = im2uint8(read(t)); % read tiff file
        im_tresh = imageData;
        im_tresh(imageData<tmpthrVal) = 1;
        im_tresh(imageData>=tmpthrVal) = 0;
        im_tresh(~tmpMask) = 0; % overlay mask on image
        img_values = im_tresh(im_tresh==1); % leave behind zero values (non-ROI pixels)
        tmpPupil(1,i) = sum(img_values); % pupil area = sum of white pixels after binarization
    end
else
    filenames = reshape(filenames,[],size(filenames,3));
    filenames(filenames<tmpthrVal) = -1;
    filenames(filenames>=tmpthrVal) = 0;
    filenames = -1*filenames;
    filenames(~tmpMask(:),:) = 0;
    tmpPupil = sum(filenames);
end

pupil_outlier = filloutliers(tmpPupil,'next','percentiles',[1 99]); % replaces outluiers with previous values
pupil_raw = 2 * ((pupil_outlier/pi).^0.5);
pupil_raw = pupil_raw./behaviorROIs.eye_length;
pupil = pupil_raw;

end