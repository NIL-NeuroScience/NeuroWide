function [tform, reg, f] = brain_register(moving, fixed)

% correct sizes
dim1 = size(fixed);
dim2 = size(moving);

if dim1 ~= dim2
    
end

% register images
[optimizer, metric] = imregconfig('multimodal');
tform = imregtform(moving, fixed, 'rigid', optimizer, metric);
reg = imwarp(moving, tform, OutputView=imref2d(dim1));

% plot registration

if nargout > 2
    f = figure;
    tiledlayout(2,2);
    
    nexttile;
    imshow(fixed / 65535);
    title('Fixed Brain');
    
    nexttile;
    imshow(moving / 65535);
    title('Moving Brain');
    
    nexttile;
    imshow(cat(3, fixed / 65535, moving / 65535, zeros(dim1)));
    title('Overlaid Brains');
    
    nexttile;
    imshow(cat(3, fixed / 65535, reg / 65535, zeros(dim1)));
    title('Overlaid Registered Brains');
end

end