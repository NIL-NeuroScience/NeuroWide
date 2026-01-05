%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               f_bpf
% author - Brad Rauscher
% 
% Performs band-pass filtering.
% 
% INPUTS:
%   sig: signal
%   fr: frequency window [f_low, f_high]. If f_low = 0, only performs
%       low-pass filtering. If f_high = fs/2, only performs high-pass 
%       filtering.
%   fs: sampling frequency of sig (Hz).
%   dim: dimension to perform filtering on. (default = 1)
% 
% OUTPUTS:
%   filt: filtered signal with the same size as sig.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function filt = f_bpf(sig,fr,fs,dim)

if nargin < 4
    dim = 1;
end

n_dims = ndims(sig);

% check dim matches number of dimensions
if dim > n_dims
    error('"sig" does not have "dim" dimensions!!!');
end

new_dims = 1:n_dims;
new_dims(new_dims == dim) = [];

new_dims = [dim, new_dims];
inv_dims = zeros(size(new_dims));
inv_dims(new_dims) = 1:numel(new_dims);

sig = permute(sig,new_dims);

% perform filtering

idx = fr == [0 fs/2];
idx = ~idx;

if idx
    [r,a] = butter(6,fr(1)/(fs/2));
    filt = filtfilt(r,a,sig);
    high = sig - filt;
    [r,a] = butter(6,fr(2)/(fs/2));
    filt = filtfilt(r,a,high);
elseif ~idx(1)
    [r,a] = butter(6,fr(2)/(fs/2));
    filt = filtfilt(r,a,sig);
elseif ~idx(2)
    [r,a] = butter(6,fr(1)/(fs/2));
    filt = filtfilt(r,a,sig);
    filt = sig - filt;
else
    filt = sig;
end

filt = permute(filt,inv_dims);

end