%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             f_coherence
% author - Brad Rauscher (created 2025)
% 
% Calculates coherence between each column in 'sig1' and 'sig2' using the
% coherencyc function in the Chronux toolbox. 
% 
% INPUTS: f_coherence(sig1, sig2, fs, _)
%   sig1: first signal (H x W x T)
%   sig2: second signal  (H x W x T)
%   fs: sampling rate of signals (Hz)
% 
% OPTIONAL INPUTS:
%   tapers: multitaper parameters [time-bandwidth product, N tapers]
%       (default = [5, 9])
% 
% OUTPUTS:
%   C: coherence value(s)
%   phi: phase value(s) for C
%   f: frequency value(s) for C
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [C, phi, f] = f_coherence(sig1, sig2, fs, varargin)
    % handle inputs
    p = inputParser;
    addParameter(p, 'tapers', [5, 9]); % multitaper parameters
    
    parse(p, varargin{:});

    % calculate coherence
    params = struct;
    params.Fs = fs;
    params.tapers = p.Results.tapers;
    params.trialave = 0;
    
    [C, phi, ~, ~, ~, f] = coherencyc(sig1, sig2, params);
    f = f';
end