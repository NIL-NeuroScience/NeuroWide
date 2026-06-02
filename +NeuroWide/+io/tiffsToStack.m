%TIFFSTOSTACK Combine all single-frame TIFFs in a folder into one multi-page TIFF stack.
%
%   tiffsToStack(inputFolder, outputFile)
%
%   inputFolder: path to the folder containing TIFF images.
%   outputFile:  path (including filename) of the output multi-page TIFF stack.
function tiffsToStack(path, outputFile)

    % Find all tiff files
    files = dir(fullfile(path, '*.tif'));
    if isempty(files)
        files = dir(fullfile(path, '*.tiff'));
    end
    if isempty(files)
        error('No TIFF files found in %s', path);
    end

    % Sort by name (alphabetical)
    [~,idx] = sort({files.name});
    files = files(idx);

    % Read first image to get size/type
    firstImg = imread(fullfile(path, files(1).name));

    % Write first image
    imwrite(firstImg, outputFile, 'tif', 'WriteMode', 'overwrite');

    % Append the rest
    for k = 2:numel(files)
        img = imread(fullfile(path, files(k).name));

        % optional check for consistency
        if ~isequal(size(img), size(firstImg))
            error('Image %s has different size than the first image', files(k).name);
        end

        imwrite(img, outputFile, 'tif', 'WriteMode', 'append');
    end

    fprintf('Saved %d frames to %s\n', numel(files), outputFile);
end