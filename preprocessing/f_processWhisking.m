function [whisker_long,whisker_pad] = f_processWhisking(filenames,behaviorROIs)

tmpwhisker = [];
tmpwhisker2 = [];

for i = 1:numel(filenames)
    t = Tiff(filenames{i},'r');
    imageData = im2uint8(read(t));
    Icropped = imcrop(imageData,behaviorROIs.whisker_rois(1,:));
    Icropped2 = imcrop(imageData,behaviorROIs.whisker_rois(2,:));
    switch i
        case 1
            img_prev = Icropped;
            img_prev2 = Icropped2;
        otherwise
            img_show = abs(Icropped - img_prev);
            img_show2 = abs(Icropped2 - img_prev2);
            if sum(img_show(:)) ~= 0
                tmpwhisker(i) = sum(img_show(:));
                tmpwhisker2(i) = sum(img_show2(:));
            end
            img_prev = Icropped;
            img_prev2 = Icropped2;
    end
end

whisker_long = rescale(tmpwhisker);
whisker_pad = rescale(tmpwhisker2);

end