function [movingPoints,fixedPoints] = f_create_landmarks(moving,fixed,N,rotate,thresh)
% Version 1.0 Harrison Fisher 12/19/2021
% create_landmarks: create pairs of points for a control point registration
%   

% plot images 
figure('Units','normalized','Position',[0.081,0.2429,0.4762,0.6085]);
ax1 = subplot(1,2,1);
imagesc(mat2gray(moving))
colormap(ax1,gray)
clim([0,thresh])
title('Alternate between picking a point on the moving image (below)')

ax2 = subplot(1,2,2);
imagesc(fixed,'AlphaData',~isnan(fixed))
colormap(ax2,parula)
title('and the corresponding point on the fixed image (below)')

% initialize points storage
movingPoints = zeros(N,2);
fixedPoints = zeros(N,2);

% alternate between moving and fixed point picking 
for n = 1:N    
    [x1,y1] = ginput(1);
    hold on
    plot(x1,y1,'r.','MarkerSize',30,'Parent',ax1)

    [x2,y2] = ginput(1);
    hold on
    plot(x2,y2,'r.','MarkerSize',30,'Parent',ax2)
    
    if rotate == 1 
        % if transposed, swap dimensions
        movingPoints(n,1) = y1;
        movingPoints(n,2) = x1;
    else
        movingPoints(n,1) = x1;
        movingPoints(n,2) = y1;      
    end
    
    fixedPoints(n,1) = x2;
    fixedPoints(n,2) = y2;
end

close all
end

