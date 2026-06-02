function [parcellation] = Finish_Registration(dat,FinalTform,tolerance,isLeft)
% Version 1.0 Harrison Fisher 12/19/2021
% Version 1.1 Patrick Doran   03/31/2022
% Load Atlas file
% Load Atlas file
atlas_file = fullfile(NeuroWide.repoPath, 'data/allen_proj.mat');
load(atlas_file);

% Load regions table
% regions_file = 'regions.csv';
% regions_table = ReformatRegionTable(regions_file);

regions_file = fullfile(NeuroWide.repoPath, 'data/regions_matlab.mat');
regions_struct = load(regions_file);
regions_table = struct2table(regions_struct);

outline_file = fullfile(NeuroWide.repoPath, 'data/allen_proj_outline.mat');
load(outline_file)

fixed = AllenAtlas;
 
moving = dat(:,:,1);

% apply transformation 
data_warped = NeuroWide.allenAtlas.ApplyRegistration(dat,fixed,FinalTform);

%% Extract parcellation timeseries 

hemi = NeuroWide.allenAtlas.Hemisphere_Mask(data_warped,isLeft);

% apply hemisphere masks to parcellation 
make_plot = 1;
parcellation = NeuroWide.allenAtlas.ApplyParcellation(data_warped,AllenAtlas,AllenOutline,regions_table,hemi,tolerance,make_plot);

% Convert Parcelation to native space
parcellation = NeuroWide.allenAtlas.TransformParcellation(parcellation,moving,FinalTform,isLeft);
% parcellation = AddRegionNames(parcellation);

%% Plot Parcels on image
% close all;
% ParcelPlot(parcellation,moving);
end