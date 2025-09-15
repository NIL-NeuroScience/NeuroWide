function metadata = f_loadAnalysis(ImagingLog,run_path,filename,varargin)

p = inputParser;
addParameter(p,'type',[]);
addParameter(p,'var',[]);
parse(p,varargin{:});

if isempty(p.Results.var)
    varName = 'metadata';
else
    varName = p.Results.var;
end

metadata = struct;
files = struct;

for Idx = 1:numel(ImagingLog)

    mouse = ImagingLog(Idx).Mouse;
    date = ImagingLog(Idx).Date;
    run = ImagingLog(Idx).Run;
    root = ImagingLog(Idx).rootFolder;
    grab = ImagingLog(Idx).GRAB;

    files.Run = fullfile(root,run_path,sprintf('Run%02i',run));
    files.load = fullfile(files.Run,filename);

    metadata_run = load(files.load,varName);
    metadata_run = metadata_run.(varName);
    
    if string(p.Results.type) == "cocaine"
        for i = 1:numel(metadata_run)
            metadata_run(i).Week = ImagingLog(Idx).Week;
            metadata_run(i).Session = ImagingLog(Idx).Session;
        end
    end
    
    dateString = strrep(date,'-','_');
    dateString = ['d' dateString];
    
    for i = 1:numel(metadata_run)
        metadata_run(i).Mouse = mouse;
        metadata_run(i).Date = date;
        metadata_run(i).Run = run;
        metadata_run(i).GRAB = grab;
    end

    metadata.(mouse).(dateString).(sprintf('Run%02i',run)) = metadata_run;
end
end