%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             f_smooth2d
% author - Brad Rauscher (created 2023)
% 
% Smooths the input frame data in 'sig' using a gaussian kernel with size
% parameter 'smooth'.
% 
% INPUTS: f_smooth2d(sig, smooth)
%   sig: input video
%   smooth: smoothing kernel size
% 
% OUTPUTS:
%   smoothed: smoothed video
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function smoothed = f_smooth2d(sig, smooth)

[nrows, ncols] = size(sig, [1, 2]);

filtim = ((1 : nrows - ceil((nrows + 1) / 2))' * ones(1, ncols)).^2 + ...
    (ones(nrows, 1) * 1 : ncols - ceil((ncols + 1) / 2)).^2;
filtim = exp(-smooth^2 * filtim / max(filtim(:)));
tmp = fftshift(fftshift(fft(fft(squeeze(sig(:, :, :)), [], 1), ...
    [], 2), 2), 1);

for fn = 1 : size(tmp, 3)
  tmp(:, :, fn) = squeeze(tmp(:, :, fn)) .* filtim;
end

smoothed = real(ifft(ifft(ifftshift(ifftshift(tmp, 2), 1), [], 2), [], 1));

end