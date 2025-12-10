function [refBM, refParcellation, refMasks] = f_loadRefAllen(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         f_loadRefAllen
% Loads reference allen parcellation variables.
% 
% OPTIONAL INPUTS:
%   hem: 0 (default) - both hemispheres. 1 - left hemisphere. 
%     2 - right hemisphere. 
%
% OUTPUTS:
%   refBM: reference brain exposure mask.
%   refParcellation: full reference parcellation struct.
%   refMasks: extracted allen region masks.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% handle inputs

p = inputParser;
addParameter(p,'hem',0);

parse(p,varargin{:});

%% load refAllen.mat

path = fullfile(f_path,'plotting/refAllen.mat');

tmp = load(path);

%% extract and process struct

refBM = tmp.refBM;
refParcellation = tmp.refParcellation;

if p.Results.hem == 2
    refParcellation.Masks = refParcellation.Masks(:,301:end,:,:);
    refBM = refBM(:,301:end,:);
elseif p.Results.hem == 1
    refParcellation.Masks = refParcellation.Masks(:,1:300,:,:);
    refBM = refBM(:,1:300,:);
end

refMasks = refParcellation.Masks;

refMasks = sum(refMasks,4);
refMasks(refMasks==0) = NaN;

end