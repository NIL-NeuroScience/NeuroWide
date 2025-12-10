function [parcellation] = f_ApplyParcellation(data_warped,AllenAtlas,AllenOutline,regions_table,hemi,tolerance,make_plot)
% Version 1.0 Harrison Fisher 12/19/2021
% ApplyParcellation: create ROI mask from the Allen atlas and extract
% timeseries from each. Each roi is first masked with the manually drawn
% hemisphere masks. The tolerance sets how much of the ROI must be retained
% after applying the hemisphere mask to include this ROI. 
%
%   Inputs:
%       data_warped:    3D array of data in template (allen atlas) space
%       AllenAtlas:     2D array of Allen Atlas      
%       regions_table:  table of Allen atlas regions 
%       hemi:           structure of hemisphere masks 
%       tolerance:      percentage of the ROI that needs to remain after
%           applying the hemisphere mask (ie .9 means at most 10% of the 
%           voxels can be lost
%       make_plot:      1 to display the image of the valid parcels, 0 to
%                       skip the plot
%   Outputs:
%       parcellation: structure containing modified ROI masks and extracted
%           timeseries. In the allTS field, NaNs are filled 

parcellation = []; 

L_count_v = 1;
L_count_inv = 1;

R_count_v = 1;
R_count_inv = 1;

T = size(data_warped,3);
N = height(regions_table);
allTS = zeros(T,2*N);

% initialize arrays for plotting
R = zeros(size(AllenAtlas));
R(AllenAtlas >= 1) = 1;
G = zeros(size(AllenAtlas));

for r = 1:N
%     disp(strcat(['Region',' ',num2str(r),'/',num2str(N)]))
    label = regions_table.acronym{r};
    
    Lparcel = poly2mask(regions_table.left_x{r}, regions_table.left_y{r},size(AllenAtlas,1),size(AllenAtlas,2));
    %Lparcel = (AllenAtlas == iAreas);
        % need to split in half if using this way??
    
    % compute overlap between manual hemisphere mask and the parcellation
    overlap = hemi.left .* Lparcel; 
    
    % check how much parcel is reduced by the hemisphere mask 
    if sum(overlap(:)) > tolerance * sum(Lparcel(:)) 
        parcellation.LeftROIs{r} = overlap;
        parcellation.LeftValid{L_count_v} = label;
        % extract timeseries
%         disp(strcat(['extracting ts for left',' ',label]))

        ts = f_ExtractTimeseries(data_warped,overlap);
        allTS(:,r) = ts;
        parcellation.LeftTS{L_count_v} = ts;
        L_count_v = L_count_v + 1;
    else
        parcellation.LeftROIs{r} = nan;
        parcellation.LeftInvalid{L_count_inv} = label;
        L_count_inv = L_count_inv + 1;
        allTS(:,r) = nan(T,1);
    end
    
    % Right Side 
    Rparcel = poly2mask(regions_table.right_x{r}, regions_table.right_y{r},size(AllenAtlas,1),size(AllenAtlas,2));
    overlap = hemi.right .* Rparcel; 
    
    if sum(overlap(:)) > tolerance * sum(Rparcel(:)) 
        parcellation.RightROIs{r} = overlap;
        parcellation.RightValid{R_count_v} = label;
%         disp(strcat(['extracting ts for right',' ',label]))

        ts = f_ExtractTimeseries(data_warped,overlap);
        parcellation.RightTS{R_count_v} = ts;

        allTS(:,N+r) = ts;
        R_count_v = R_count_v + 1;
    else
        parcellation.RightROIs{r} = nan;
        parcellation.RightInvalid{R_count_inv} = label;
        R_count_inv = R_count_inv + 1;
        
        allTS(:,N+r) = nan(T,1);
        
    end
    
    % make array for plotting 
    R(parcellation.LeftROIs{r} == 1) = 0;
    G(parcellation.LeftROIs{r} == 1) = 1;
    
    R(parcellation.RightROIs{r} == 1) = 0;
    G(parcellation.RightROIs{r} == 1) = 1; 
end

parcellation.allTS = allTS; 

brainmask = (AllenAtlas >= 1);

regions_img = zeros([size(AllenAtlas),3]);
regions_img(:,:,1) = R;
regions_img(:,:,2) = G;

if make_plot == 0
    imagesc(regions_img,'AlphaData',brainmask)
    plot_allen_outlines(regions_table,AllenOutline,'black')
end

parcellation.regions_img = regions_img;
parcellation.brainmask = brainmask;

end
