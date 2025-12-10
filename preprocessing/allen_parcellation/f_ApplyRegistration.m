function [out] = f_ApplyRegistration(data,template,tform)
% Version 1.0 Harrison Fisher 12/19/2021
% Version 1.1 Patrick Doran   03/31/2022
%ApplyRegistration: apply the transformation matrix to every volume of the
% 3D image
%   Inputs:
%       data: Npixels x Mpixels x timepoints 3D matrix of data
%       template: 2D array of target standard space atlas 
%       tform: affine or similarity transformation matrix to apply to each volume 
%
%
%
N = size(data,3);
out_dims = size(template);
ref = imref2d(out_dims);

out = zeros([out_dims,N]);

for t = 1:N
    out(:,:,t) = imwarp(data(:,:,t),tform,'OutputView',ref);
    
end

end