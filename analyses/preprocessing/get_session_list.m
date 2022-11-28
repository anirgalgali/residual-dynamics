function [idlist_valid, idlist_invalid] = get_session_list(data_dir, animal, varargin)
%{ Extracts the list of recording sessions from a list of data files, where 
%  each data file is a single recording session
% Input
% - data_dir(struct) - contianing information about the directory where the
%   data is stored
% - animal (string) -identifier associated with the dataset, typically the
%   name of the animal
% -varargin
%  - sessions_to_exclude (cell array) - each element of the cell array is a
%  string label indicating which files to exclude while populating the
%  session list.
% Output
% - idlist_valid (cell array) - provides the name of all the "valid"
%   sessions"
% - idlist_invalid (cell array) - all the excluded sessions
% Author : Aniruddh Galgali (Jan 2022)
%}
if(nargin == 3)
    sessions_to_exclude = varargin{1};
elseif(nargin == 2)
    sessions_to_exclude = {};
end

idx_animal = false(1,length(data_dir.idlist));
for iexp = 1:length(data_dir.idlist)   
    idx_animal(iexp) = any(strfind(data_dir.idlist{iexp},[animal]));   
end
idlist_animal = data_dir.idlist(idx_animal);

n_sessions = length(idlist_animal);

if(length(sessions_to_exclude) >= 1)
    invalid_session_indices = NaN(length(sessions_to_exclude),1);
    for iexp = 1: length(sessions_to_exclude)    
        invalid_session_indices(iexp) = find(cell2mat(cellfun(@(x) ~isempty(x),strfind(idlist_animal, sessions_to_exclude{iexp}),'uni',false)));
    end
    idlist_valid = idlist_animal(setdiff(1:n_sessions,invalid_session_indices));
    idlist_invalid = idlist_animal(invalid_session_indices);
else
    idlist_valid = idlist_animal;
    idlist_invalid = [];
end

end