function f_savePNG(filename,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         f_savePNG
% author - Brad Rauscher (created Dec/18/2025)
% 
% Saves current figure as .png
% 
% INPUTS:
%   filename: name to save file as.
%   flag: (optional) set to 0 to not save.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% handle inputs

p = inputParser;
addParameter(p,'dim',[]);
addParameter(p,'flag',1);

parse(p,varargin{:});

%%

ppi = get(groot,"ScreenPixelsPerInch");

if p.Results.flag
    if isempty(p.Results.dim)
        exportgraphics(gcf,filename,Resolution=300,BackgroundColor='white');
    else
        exportgraphics(gcf,filename,Resolution=300,BackgroundColor='white', Units='inches', Width=p.Results.dim(2), Height=p.Results.dim(1));
    end
end

end