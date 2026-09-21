function reg_files = find_reg(file_name, var_name, root, level)

% gather all files
files = {root};

for i = 1:level
    new_files = {};
    for f = 1:numel(files)
        tmp_files = dir(files{f});
        tmp_files(1:2) = [];
        new_files = [new_files; fullfile(files{f}, {tmp_files.name}')];
    end
    files = new_files;
end
N = numel(files);

% check for file_name
file_exists = false(N, 1);
for i = 1:N
    file_exists(i) = logical(exist(fullfile(files{i}, file_name)));

    if file_exists(i) & ~isempty(var_name)
        var_info = whos('-file', fullfile(files{i}, file_name));
        if ~sum(strcmp(var_name, {var_info.name}))
            file_exists(i) = false;
        end
    end
end

reg_files = fullfile(files(file_exists), file_name);

end