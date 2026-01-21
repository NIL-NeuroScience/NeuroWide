%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              f_spectra
% author - Brad Rauscher (created 2025)
% 
% Calculates power spectral density for each column in 'sig' using the
% mtspectrumc function in the Chronux toolbox.
% 
% INPUTS: f_spectra(sig, fs, varargin)
%   sig: signal
%   fs: sampling rate of signals (Hz)
% 
% OPTIONAL INPUTS:
%   tapers: multitaper parameters [time-bandwidth product, N tapers]
%       (default = [5, 9])
% 
% OUTPUTS:
%   S: spectral power density value(s)
%   f: frequency value(s) for S
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [S, f] = f_spectra(sig, fs, varargin)
    % handle inputs
    p = inputParser;
    addParameter(p, 'tapers', [5, 9]); % multitaper parameters
    
    parse(p, varargin{:});

    % calculate coherence
    params = struct;
    params.Fs = fs;
    params.tapers = p.Results.tapers;
    params.trialave = 0;
    
    [S, f] = mtspectrumc(sig, params);
    f = f';
end