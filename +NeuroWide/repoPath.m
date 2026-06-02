%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              repoPath
% author - Brad Rauscher (created 2026)
% 
% Returns parent directory of function.
% 
% OUTPUTS:
%   path: parent directory of function
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function path = repoPath()

path = fileparts(fileparts(mfilename('fullpath')));

end