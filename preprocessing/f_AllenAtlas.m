function [parcellation] = f_AllenAtlas(img,brain_mask,isLeft)
%%

Tolerance = 0.01;

tmpIn = img;
% tmpIn = uint16(tmpIn);

f_Do_Registration(tmpIn,brain_mask);

tmpIn = double(tmpIn).*brain_mask;
% tmpIn = uint16(tmpIn);

parc = f_Finish_Registration(tmpIn,evalin('base','FinalTform'),Tolerance,isLeft);

%% find indexes of chosen areas
areas = {'MOp','MOs','SSp-bfd','SSp-tr','SSp-ll','SSp-ul','SSp-un','VISpm','VISrl','VISam','VISa','VISp'};
a = struct;

for i = 1:size(areas,2)
    try
        a.indL(i) = find(strcmp(parc.LeftValid,areas(i)));
    catch
        a.indL(i) = 0;
    end
end

for i = 1:size(areas,2)
    try
        a.indR(i) = find(strcmp(parc.RightValid,areas(i)));
    catch
        a.indR(i) = 0;
    end
end

% if sum(a.indL==0) > 0
%     error('Redo registration. Does not contain all areas!')
% end
% if sum(a.indR==0) > 0
%     error('Redo registration. Does not contain all areas!')
% end

%%

parcellation.labelsLR = areas;

if isLeft == 2
    parcellation.Masks = parc.LeftROIs(:,:,a.indL);
    parcellation.Masks(:,:,:,2) = parc.RightROIs(:,:,a.indR);

elseif isLeft == 1
    parcellation.Masks = parc.LeftROIs(:,:,a.indL);
    parcellation.Masks(:,:,:,2) = zeros(size(parcellation.Masks));

elseif isLeft == 0
    parcellation.Masks = zeros(size(parc.RightROIs,1),size(parc.RightROIs,2),size(areas,2));
    parcellation.Masks(:,:,:,2) = parc.RightROIs(:,:,a.indR);
    
end


end