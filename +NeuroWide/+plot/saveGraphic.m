function saveGraphic(path, varargin)

p = inputParser;
addParameter(p, 'save_dir', '/projectnb/devorlab/bcraus');
addParameter(p, 'flag', 1);
addParameter(p, 'dim', []);
addParameter(p, 'as_svg', 0);
parse(p, varargin{:});

fullpath = fullfile(p.Results.save_dir, path);

[dir,name,~] = fileparts(fullpath);

[~,~,~] = mkdir(dir);

filename = fullfile(dir, [name,'.png']);

if p.Results.flag
    if isempty(p.Results.dim)
        exportgraphics(gcf, filename, ...
            Resolution=300, ...
            BackgroundColor='white');
    else
        exportgraphics(gcf, filename, ...
            Resolution=300, ...
            BackgroundColor='white', ...
            Units='inches', ...
            Width=p.Results.dim(2), ...
            Height=p.Results.dim(1));
    end
end

if p.Results.as_svg
    print(gcf, fullfile(dir, [name,'.svg']), '-dsvg');
end

end