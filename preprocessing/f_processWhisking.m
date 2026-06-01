function [whisker_long,whisker_pad] = f_processWhisking(filenames,behaviorROIs)

tmpwhisker = [];
tmpwhisker2 = [];

if iscell(filenames)
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
else
    b = behaviorROIs.whisker_rois;
    b(:,3:4) = ceil(b(:,1:2) + b(:,3:4));
    b(:,1:2) = floor(b(:,1:2));

    img1 = filenames(b(1,2):b(1,4),b(1,1):b(1,3),:);
    img2 = filenames(b(2,2):b(2,4),b(2,1):b(2,3),:);
    
    tmpwhisker = std(img1(:,:,2:end) - img1(:,:,1:end-1),0,[1,2]);
    tmpwhisker2 = std(img2(:,:,2:end) - img2(:,:,1:end-1),0,[1,2]);
    tmpwhisker = [0,tmpwhisker(:)'];
    tmpwhisker2 = [0,tmpwhisker2(:)'];

end
whisker_long = rescale(tmpwhisker);
whisker_pad = rescale(tmpwhisker2);

end