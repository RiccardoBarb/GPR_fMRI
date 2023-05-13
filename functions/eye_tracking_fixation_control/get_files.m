function [folders] = get_files(directory, search_index)

% this function is written to get a list of folders from a certain
% directory. Plus it will look in the folders if there are *.asc files.
% These files are created from *.edf files that contain our eye tracking
% data. The current specific indecees for the *.asc files is L1 (localizer first load)
% and T1 (task first load) the function will also give mark if both files
% are existing

% in:   - directory:        directory where the folder are located
%       - search_index:     common string in the folders to separate them from others
% out:  - folders:          structure containing name of folders, number of
%                           *.asc files in that folder and number of
%                           *L1.asc and *T1.asc files in that folder

    
    content = dir(directory);
    name = {content.name};
    list = contains(name, search_index); % because all subject id's start with s3xx
    subj_folders=(name(list==1));

    field='subject';value= subj_folders;
    folders=struct(field,value);



    for n=1:length(folders)
        s_dir=fullfile(directory,folders(n).subject);
        cont_dir=[dir(s_dir)];
        cont_dir={cont_dir.name};
        folders(n).is_asc=sum(contains(cont_dir,'asc'));
        
        %is there a task asc file and if yes save a directory
        folders(n).is_Tasc=sum(contains(cont_dir,'T1.asc'));
        if folders(n).is_Tasc==1
            file=[folders(n).subject(1:end-2) 'd' folders(n).subject(end) 'T1.asc'];
            folders(n).Tasc_dir=fullfile(directory, folders(n).subject, file);
        else
            folders(n).Tasc_dir=NaN;
            fprintf('Error: %s does not contain a task eyetracking file\n',folders(n).subject)
        end
        
        folders(n).is_mat=sum(contains(cont_dir,'mri'));
        if folders(n).is_mat==1
            file=[folders(n).subject(1:end-2) '_mri_day' folders(n).subject(end) '.mat'];
            folders(n).mat_dir=fullfile(directory, folders(n).subject, file);
        else
            folders(n).mat_dir=NaN;
            fprintf('Error: %s does not contain an task mat file\n',folders(n).subject)
        end
    end

