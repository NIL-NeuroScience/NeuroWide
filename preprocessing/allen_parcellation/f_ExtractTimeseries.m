function [ts] = f_ExtractTimeseries(data_warped,mask)
% ExtractTimeseries: take average over the mask region for each timepoint
%   Inputs:
%       data_warped:    input 3D data 
%       mask:           binary mask (2D) of region of interest

T = size(data_warped,3);
ts = zeros([T,1]);

for t = 1:T
   tmp = data_warped(:,:,t);
   ts(t) = mean(tmp(mask==1));
end

end

