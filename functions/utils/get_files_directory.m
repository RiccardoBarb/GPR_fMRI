function [files_dir, subj_name]=get_files_directory(folder_dir,search_str,file_extension)


% this function is written to search for files in a certain directory and
% output their path.

% IN:
% folder_dir:       strings in cell structure eg {'path_a/esperiment1'},{'path_b/experiment'},.. also a
%                   single folder dir must be written in a cell.
% search_str:       string, eg name of the file eg 's01'
% file_extension:   string, file extension of the file eg '.mat'

% if size(folder_dir,1)==1
%     folder_dir=folder_dir{1};
% end
    

files_dir=cell(size(folder_dir,1),1);
subj_name=cell(length(folder_dir),1);

for fold_=1:size(folder_dir,1)

    %get directory
    files=dir(cell2mat(folder_dir(fold_)));
    name={files.name};
    search_list = contains(name, search_str );
    extent_list = contains(name, file_extension );
    list=search_list & extent_list;
    files_list=find(list==1);
    
%    subj_name(fold_)={files(files_list(1)).name(1:4)};
    
    for file_=1:length(files_list)
        
    files_dir(fold_,file_)={fullfile(files(files_list(file_)).folder,files(files_list(file_)).name)};
    subj_name(fold_,file_)={files(files_list(file_)).name(1:4)};
    end
    % get name
%    subj_name(fold_)={files(list).name(1:4)};
    
end
   
   
empt=cellfun(@isempty,files_dir);
files_dir=files_dir(~empt);
subj_name=unique(subj_name(~empt));

fprintf('\n found %i files:\n',length(files_dir))
disp(files_dir)
