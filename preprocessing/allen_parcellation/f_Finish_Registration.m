function [parcellation] = f_Finish_Registration(dat,FinalTform,tolerance,isLeft)
% Version 1.0 Harrison Fisher 12/19/2021
% Version 1.1 Patrick Doran   03/31/2022
% Load Atlas file
% Load Atlas file
atlas_file = 'allen_proj.mat';
load(atlas_file);

% Load regions table
% regions_file = 'regions.csv';
% regions_table = f_ReformatRegionTable(regions_file);

regions_file = 'regions_matlab.mat';
regions_struct = load(regions_file);
regions_table = struct2table(regions_struct);

outline_file = 'allen_proj_outline.mat';
load(outline_file)

fixed = AllenAtlas;
 
moving = dat(:,:,1);

% apply transformation 
data_warped = f_ApplyRegistration(dat,fixed,FinalTform);

%% Extract parcellation timeseries 

hemi = f_Hemisphere_Mask(data_warped,isLeft);

% apply hemisphere masks to parcellation 
make_plot = 1;
parcellation = f_ApplyParcellation(data_warped,AllenAtlas,AllenOutline,regions_table,hemi,tolerance,make_plot);

% Convert Parcelation to native space
parcellation = f_TransformParcellation(parcellation,moving,FinalTform,isLeft);
% parcellation = f_AddRegionNames(parcellation);

%% Plot Parcels on image
% close all;
% f_ParcelPlot(parcellation,moving);
end