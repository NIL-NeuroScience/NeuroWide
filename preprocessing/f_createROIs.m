function [brain_mask,hem] = f_createROIs(template)

templateLim = prctile(template(:),[1 99]);

rFig = figure;
imagesc(template);
axis image off;
colormap cmpinf;
clim(templateLim);

hem = input('Number of hemispheres? (2-both, 1-left, 0-right)');

title('Draw 1st Hemisphere');
hem1 = drawpolygon;
hem1 = createMask(hem1);
hem1 = double(hem1);hem1(hem1 == 0) = NaN;
if hem > 1
    imagesc(template);
    axis image off;
    clim(templateLim);
    title('Draw 2nd Hemisphere');
    hem2 = drawpolygon;
    hem2 = createMask(hem2);
    hem2 = double(hem2);hem2(hem2 == 0) = NaN;
    brain_mask = hem2;
    brain_mask(~isnan(hem1)) = 1;
else
    brain_mask = hem1;
end
close(rFig);

end