function [Image_Parcellation] = f_TransformParcellation(parcellation,image,transform,isLeft)
% Version 1.0 Harrison Fisher 12/19/2021
% TransformParcellation transforms parcellation from atlas space to image
% space
%   Inputs:
%    Parcellation: A structure that includes the parcellation in atlas
%    space
%    image: Any image in the image space dimensions that we need to
%    transform to
%    transform: The affine transformation that was used to transform the
%    image to atlas space
%   Outputs:
%      Image_Parcellation: The parcellation in image space with ROIs that
%      define which image pixels are in which atlas region

Image_Parcellation = [];

if isLeft == 2
    Image_Parcellation.LeftValid = parcellation.LeftValid;
    Image_Parcellation.RightValid = parcellation.RightValid;
    LeftROIs = {};
    RightROIs = {};
    
    counter = 1;
    for x = 1:length(parcellation.LeftROIs)
        if ~isnan(cell2mat(parcellation.LeftROIs(x)))
            LeftROIs(counter) = parcellation.LeftROIs(x);
            counter = counter + 1;
        end
    end
    
    counter = 1;
    for x = 1:length(parcellation.RightROIs)
        if ~isnan(cell2mat(parcellation.RightROIs(x)))
            RightROIs(counter) = parcellation.RightROIs(x);
            counter = counter + 1;
        end
    end
    
    LeftImageROIs = zeros(size(image,1),size(image,2),length(LeftROIs));
    RightImageROIs = zeros(size(image,1),size(image,2),length(RightROIs));
    
    InverseTransform = invert(transform);
    
    for x = 1:length(LeftROIs)
        temp = f_ApplyRegistration(cell2mat(LeftROIs(x)),image,InverseTransform);
        temp(find(temp>=0.5)) = 1;
        temp(find(temp<0.5)) = 0;
        LeftImageROIs(:,:,x) = temp;  
    end
    
    LeftImageROIs = logical(LeftImageROIs);
    
    for x = 1:length(RightROIs)
        temp =f_ApplyRegistration(cell2mat(RightROIs(x)),image,InverseTransform);
        temp(find(temp>=0.5)) = 1;
        temp(find(temp<0.5)) = 0;
        RightImageROIs(:,:,x) = temp;
    end
    
    RightImageROIs = logical(RightImageROIs);
    
    Image_Parcellation.LeftROIs = LeftImageROIs;
    Image_Parcellation.RightROIs = RightImageROIs;

elseif isLeft == 1
    Image_Parcellation.LeftValid = parcellation.LeftValid;
    LeftROIs = {};
    
    counter = 1;
    for x = 1:length(parcellation.LeftROIs)
        if ~isnan(cell2mat(parcellation.LeftROIs(x)))
            LeftROIs(counter) = parcellation.LeftROIs(x);
            counter = counter + 1;
        end
    end
    
    LeftImageROIs = zeros(size(image,1),size(image,2),length(LeftROIs));
    
    InverseTransform = invert(transform);
    
    for x = 1:length(LeftROIs)
        temp = f_ApplyRegistration(cell2mat(LeftROIs(x)),image,InverseTransform);
        temp(find(temp>=0.5)) = 1;
        temp(find(temp<0.5)) = 0;
        LeftImageROIs(:,:,x) = temp;  
    end
    
    LeftImageROIs = logical(LeftImageROIs);
    
    Image_Parcellation.LeftROIs = LeftImageROIs;
elseif isLeft == 0
    Image_Parcellation.RightValid = parcellation.RightValid;
    RightROIs = {};
    
    counter = 1;
    for x = 1:length(parcellation.RightROIs)
        if ~isnan(cell2mat(parcellation.RightROIs(x)))
            RightROIs(counter) = parcellation.RightROIs(x);
            counter = counter + 1;
        end
    end
    
    RightImageROIs = zeros(size(image,1),size(image,2),length(RightROIs));
    
    InverseTransform = invert(transform);
    
    for x = 1:length(RightROIs)
        temp =f_ApplyRegistration(cell2mat(RightROIs(x)),image,InverseTransform);
        temp(find(temp>=0.5)) = 1;
        temp(find(temp<0.5)) = 0;
        RightImageROIs(:,:,x) = temp;
    end
    
    RightImageROIs = logical(RightImageROIs);
    
    Image_Parcellation.RightROIs = RightImageROIs;
end

end