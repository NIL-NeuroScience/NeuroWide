function parcellation = parcellate(sig, masks, varargin)

p = inputParser;
addParameter(p, 'type', 'mean');
parse(p, varargin{:});

masks = reshape(masks, [], size(masks, 3));

masks(isnan(masks)) = 0;

dim = size(sig);

sig = reshape(sig, [], dim(3))';
sig(isnan(sig)) = 0;

parcellation = sig * masks;
if strcmp(p.Results.type, 'mean')
    parcellation = parcellation ./ sum(masks, 1);
end

end