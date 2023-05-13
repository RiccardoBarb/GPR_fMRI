function [folder_dir]=get_data_folder_directory(data_dir,search_index, exclude)%,varargin)

% this function searches in a directory for folders with the given search
% index.This is very practical if you search for certain folders that
% contain data that you want to process, while other folders you do want to
% leave asside. eg process subject a and b, but leave subject c and d
% asside.
%
% IN:   data_dir:       string, directory to search in
%       search_index:   string, which string to search for in the folder
%                       names
%       exclude:        optional, Cell array, each cell containing a string that are
%                       not of interest eg: exclude={'s310','s311'};
%                       
%
% Out:  foder_dir:      cell array, each cell contains the fullpath of the
%                       folder that contains the search_index found in the
%                       given directory

if nargin<3
    exclude=[];
end


%%%%%%%%
% get list of subjects
folders=dir(data_dir);
name={folders.name};
%list = contains(name, search_index );

    list=strfind(name,search_index);
    listn=zeros(size(list));
    for n=1:length(list)
        if isempty(list{n})
            listn(n)=0;
        else
            if isfolder(fullfile(folders(n).folder,folders(n).name))
                listn(n)=1;
            else
                listn(n)=0;
                name{n}='.';
            end
        end
    end

subj_list=(name(listn==1));
fprintf('\n%i folders found\n\n',length(subj_list))

% exclude subjects

if ~isempty(exclude)
    for n=1:length(exclude)
        index = find(contains(subj_list,exclude(:)));
        

    end
        subj_list(index)=[];
else
    index=[];
end
fprintf('\n%i folders excluded\n\n',length(index))
fprintf('\n%i folders will be processed\n\n',length(subj_list))

% get directory
% allocate
folder_dir=cell(length(subj_list),1);



for f_=1:length(subj_list)
    
    subj=(strfind(name,subj_list(f_)));
    id=find(~cellfun(@isempty,subj));
    folder_dir(f_)={fullfile(folders(id).folder,folders(id).name)};
    
end