%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            f_downsample
% author - Brad Rauscher (created 2024)
% 
% Downsamples the frame data in 'sig' into 'bin' sized square bins.
% 
% INPUTS: f_downsample(sig, bin)
%   sig: signal to downsample
%   ds: downsampling factor, size of square bin
% 
% OUTPUTS:
%   downsamples: downsampled signal
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function downsampled = f_downsample(sig, ds)

if ds == 1
    downsampled = sig;
else   
    dim = size(sig);
    if size(dim, 2) == 2
        dim(3) = 1;
    end
    downsampled = zeros(floor(dim(1) / ds), floor(dim(2) / ds), dim(3));
    
    for h = ds : ds : dim(1)
        for w = ds : ds : dim(2)
            downsampled(h / ds, w / ds, :) = ...
                mean(sig(h - ds + 1 : h, w - ds + 1 : w, :), [1, 2]);
        end
    end
end
