%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              f_path
% author - Brad Rauscher (created 2025)
% 
% Returns parent directory of function.
% 
% OUTPUTS:
%   path: parent directory of function
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function path = f_path()

path = fileparts(mfilename('fullpath'));

end