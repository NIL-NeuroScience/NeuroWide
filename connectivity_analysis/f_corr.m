function [corr,f] = f_corr(sig1,sig2,dim,varargin)

% parse inputs
p = inputParser;
addParameter(p,'plot',false);

parse(p,varargin{:});

% calculate correlation along dim
sig1 = sig1 - mean(sig1,dim);
sig2 = sig2 - mean(sig2,dim);

std1 = sum(sig1.^2,dim);
std2 = sum(sig2.^2,dim);
std = sqrt(std1.*std2);
cov = sum((sig1.*sig2),dim);

corr = cov./std; 

% plot result
if p.Results.plot
    f = figure;
    imagesc(corr,AlphaData=~isnan(corr));
    axis image off;
    c = colorbar;
    c.Label.String = 'r';
    set(gca,FontSize=14);
else
    f = [];
end