%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             f_2Dconvmtx
% author - Brad Rauscher (created Jan 16th, 2026)
% 
% Extends the functionality of MATLAB's convmtx function to convolve each
% column of the 2D matrix h. 
% 
% INPUTS: f_2Dconvmtx(h, n)
%   h: T x N matrix to convolve
%   n: n length vector to convolve with
% 
% OUTPUTS:
%   A: convolution product
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function A = f_2Dconvmtx(h, n)

N = size(h, 2);

additional_zeros = zeros(n, N);

m = size(h, 1) + n - 1;
A = repmat([h; additional_zeros], 1, 1, n + 1);

A = permute(A, [1, 3, 2]);
A = reshape(A, [], N);
A = reshape(A(1 : m * n, :), m, n, N);

end