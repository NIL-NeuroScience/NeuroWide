function [mat] = f_nan2zero(mat)

mat(isnan(mat)) = 0;
mat = logical(mat);