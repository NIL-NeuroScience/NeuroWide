function [] = f_Do_Registration(dat,brain_mask)

% Load Atlas file
script_path = fileparts(mfilename('fullpath'));

atlas_file = fullfile(script_path,'allen_proj.mat');
load(atlas_file);

% Load regions table
regions_file = fullfile(script_path,'regions_matlab.mat');
regions_struct = load(regions_file);
regions_table = struct2table(regions_struct);

outline_file = fullfile(script_path,'allen_proj_outline.mat');
load(outline_file)

fixed = AllenAtlas;

time_series2 = dat.*brain_mask;
rangeDat = prctile(time_series2(:),[1 99]);
moving = 65535 *(dat-rangeDat(1))/range(rangeDat);
moving(moving<0) = 0;
moving(moving>65535) = 65535;

% error('tmp');
%% Registration 
% set registration parameters
    % similarity - rotation, scaling, and translation
    % affine     - similarity + shearing 
RegType = 'similarity';  
rotate_first = 0;   % set to 1 to automatically rotate 90 degrees first 
N_points = 2;       % number of points to use in landmark registration 
threshold = 1;     % grayscale threshold to increase contrast for registration  
New_Registration = 1; % Turns to 0 when user is satisfied with registration


    % pick landmark / control points to speed up the manual registration
    InitialTform = f_InitialRegistration(moving,fixed,N_points,rotate_first,'nonreflectivesimilarity',threshold);
        % this function computes an initial transformation based on the control
        % points that can be editted with the app below:
%     moving = moving.*brain_mask;
    % Refine registration interactively 
    app = InteractiveImageRegistration(moving,fixed,regions_table,AllenOutline,InitialTform,threshold);

    while isvalid(app)
        pause(0.1);
    end

    
end
