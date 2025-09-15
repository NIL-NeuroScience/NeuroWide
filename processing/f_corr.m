function [corr,f] = f_corr(varargin)

% parse inputs
sig1 = varargin{1};
sig2 = varargin{2};
dim = varargin{3};

p = inputParser;
addParameter(p,'plot',false);

parse(p,varargin{4:end});

sig1 = sig1 - mean(sig1,dim);
sig2 = sig2 - mean(sig2,dim);

std1 = sum(sig1.^2,dim);
std2 = sum(sig2.^2,dim);
std = sqrt(std1.*std2);
cov = sum((sig1.*sig2),dim);

corr = cov./std; 

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