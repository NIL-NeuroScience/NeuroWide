function ROIs = drawBehaviorROIs(imageData)

bFig = figure;
imagesc(imageData);
clim([0, 255]);
colormap gray;
axis image off;

title('Draw ellipse around the eye');
roi = drawellipse('Color','r'); % for best results put start and end points in the corners of the eye 
pause;
ROIs.pupil_mask = createMask(roi);
delete(roi);

title('Draw point in the pupil');
pupilpos = drawpoint('Color','r'); % for best results put start and end points in the corners of the eye 
pause;
point.pupil.pos = pupilpos.Position;
pval = imageData(round(pupilpos.Position(2)),round(pupilpos.Position(1)));
delete(pupilpos);

title('Draw point in the iris');
irispos = drawpoint('Color','r'); % for best results put start and end points in the corners of the eye 
pause;
point.iris.pos = irispos.Position;
ival = imageData(round(irispos.Position(2)),round(irispos.Position(1)));
ROIs.pupil_thrVal = (pval+ival)/2;
delete(irispos);

title('Draw line across eye (horizontal)');
eyelength = drawline('Color','r');
pause;
ROIs.eye_length = ((eyelength.Position(1,1)-eyelength.Position(2,1))^2+(eyelength.Position(1,2)-eyelength.Position(2,2))^2)^0.5;
delete(eyelength);

title('Select ROI around long whiskers');
r = drawrectangle('Color','r');
roi1 = r.Position;
pause;
delete(r);

title('Select ROI around whisker pad');
r2 = drawrectangle('Color','r');
roi2 = r2.Position;
pause;
ROIs.whisker_rois = [roi1; roi2];

close(bFig);

end