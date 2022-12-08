function project = ProjectLoad_v2(projectname);
%
% Load an analysis project file 
%
% project = ProjectLoad(projectname);
%

global DIRS
project = [];

projectpath = fullfile(DIRS.analysis,projectname);

if exist(projectpath)~=7
    disp('project does not exist');
    return;
end

projectfn = fullfile(projectpath,[projectname '.mat']);

if exist(projectfn)==2
    load(projectfn);
else
    disp('could not load project file. does not exist.');
end

% if ~isfield(project,'idlist')
%     project.idlist = ProjectCreateIdList(projectname);
% end
project.idlist = ProjectCreateIdList(projectname);

return

function newidlist = ProjectCreateIdList(projectname)
global DIRS

    filelist = dir(fullfile(DIRS.analysis,projectname));
    fnlist = {filelist.name};
    fnlist = setdiff(fnlist,{'.','..','tools','metadata','.DS_Store',[projectname '.mat']});
    
    ntrim = length('.mat');

    newidlist={};
    for ilist = 1:length(fnlist)
        str = fnlist{ilist}(1:end-ntrim);
        % str(findstr(str,'_'))='.';        
        newidlist{ilist} = str;
    end

return