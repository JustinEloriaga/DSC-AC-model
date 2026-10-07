function directory=dsc_create_timestamped_directory(parent,prefix,suffix)
% Create a second-resolution timestamp directory without a numeric suffix.
if ~exist(parent,'dir'), mkdir(parent); end
while true
    stamp=char(datetime('now','Format','yyyyMMdd-HHmmss'));
    directory=fullfile(parent,[char(prefix) stamp char(suffix)]);
    file=java.io.File(directory);
    if file.mkdir()
        return
    end
    if ~isfolder(directory)
        error('dsc:RunDirectory','Could not create run directory: %s',directory);
    end
    pause(0.05);
end
end
