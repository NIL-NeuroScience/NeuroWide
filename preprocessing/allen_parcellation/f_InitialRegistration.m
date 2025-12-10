function [initial_tform] = f_InitialRegistration(moving,allen,Ncp,rotate,transform_type,threshold)
% Version 1.0 Harrison Fisher 12/19/2021
% Version 1.1 Sreekanth Kura  07/19/2023
% InitialRegistration Perform an initial registration of data to the allen
% atlas 
%   Inputs:
%       moving:         Actual data collected from imaging system
%       allen:          Allen atlas template
%       Ncp:            number of control points to use
%       rotate:         1 to rotate image 90 initially, 0 to leave as is
%       transform_type: method of transformation to be used by fitgeotrans
%                        ie 'similarity','affine' etc. 
if rotate == 1
    moving_t = moving';
else 
    moving_t = moving;
end 


[movingPoints,~] = f_create_landmarks(moving_t,allen,Ncp,rotate,threshold);

% Here we define bregma and lambda as the fixed points. This means the
% points the user chooses on the Allen Atlas will be disregarded.
fixedPoints =[115.1806  103.8459;
  115.1806  196.5870];

landmark_tform = fitgeotrans(movingPoints,fixedPoints,transform_type);

initial_tform = landmark_tform;
end

