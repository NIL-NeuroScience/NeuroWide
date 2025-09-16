function [pupil,whisker_pad,whisker_long] = f_calcBehavior(behCam,behaviorROIs)

tmpMask = behaviorROIs.pupil_mask;
tmpthrVal = behaviorROIs.pupil_thrVal;

pupilROI = reshape(behCam,[],size(behCam,3));
pupilROI = pupilROI(tmpMask,:);
pupilROI(pupilROI<tmpthrVal) = 0;
pupilROI(pupilROI>=tmpthrVal) = 1;
pupilROI = sum(pupilROI==0,1);

pupil = filloutliers(pupilROI,'next','percentiles',[0.5 99.5]); % replaces outluiers with previous values
pupil = 2*((pupil/pi).^0.5);
pupil = pupil./behaviorROIs.eye_length;

behaviorROIs.whisker_rois = round(behaviorROIs.whisker_rois);

ROI_long = behCam(behaviorROIs.whisker_rois(1,2):behaviorROIs.whisker_rois(1,2)+behaviorROIs.whisker_rois(1,4),behaviorROIs.whisker_rois(1,1):behaviorROIs.whisker_rois(1,1)+behaviorROIs.whisker_rois(1,3),:);
ROI_pad = behCam(behaviorROIs.whisker_rois(2,2):behaviorROIs.whisker_rois(2,2)+behaviorROIs.whisker_rois(2,4),behaviorROIs.whisker_rois(2,1):behaviorROIs.whisker_rois(2,1)+behaviorROIs.whisker_rois(2,3),:);

ROI_long = ROI_long(:,:,1:end-1)-ROI_long(:,:,2:end);
ROI_pad = ROI_pad(:,:,1:end-1)-ROI_pad(:,:,2:end);

N = size(ROI_long,3);

whisker_long = rescale([0 std(reshape(ROI_long,[],N),0,1)]');
whisker_pad = rescale([0 std(reshape(ROI_pad,[],N),0,1)]');

end