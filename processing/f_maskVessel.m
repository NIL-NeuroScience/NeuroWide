function vesselMask = f_maskVessel(img,sensitivity,brain_mask,maxLength)

imAlpha = brain_mask;
imAlpha(isnan(brain_mask)) = 0;
img = img.*imAlpha;

J = fibermetric(img,sensitivity,'ObjectPolarity','dark');

% Apply segmentation. Maybe thersholding might work best
vesselMask = imbinarize(J,'adaptive');

% Do connectivity analysis and remove noisy and unconnected segements
CC = bwconncomp(vesselMask);

if nargin < 4
    maxLength = 6000;
end

for u = 1:length(CC.PixelIdxList)
    if length(CC.PixelIdxList{u}) <=300 || length(CC.PixelIdxList{u}) >= maxLength
        vesselMask(CC.PixelIdxList{u}) = 0;
    end
end

vesselMask = ~vesselMask;
vesselMask = double(vesselMask);
vesselMask(vesselMask==0) = NaN;
vesselMask = vesselMask.*brain_mask;