function bool = AnalysisDelete(projectname,idlist);
% AnalysisDelete - Erase existing analysis from a project
%
% function bool = AnalysisDelete(projectname,unitid);
%

global DIRS;
project = ProjectLoad(projectname);
projectpath = fullfile(DIRS.analysis,project.name);

bool = 0;

if ~iscell(idlist);idlist={idlist};end;

if ~strcmp(project.author,DIRS.user)
    warning('you do not own this project. analysis delete aborted.');
    return
end

if length(idlist)>1
    disp(idlist);
    prompt = sprintf('Are you sure you want to delete %s for these units (yes or no)?',project.name);
    str = input(prompt,'s');
    if strcmpi(str,'no');
        return;
        disp('aborted');
    end    
end

for iunit = 1:length(idlist)
    id=idlist{iunit};
	analysisfn = [analysisfullpathname(projectname,id) '.mat'];
	if exist(analysisfn)
		delete(analysisfn);
        fprintf(1,'%s %s deleted \n',projectname,id);
		bool = 1;
	else
        fprintf(1,'%s no analysis found\n',id);
	end
end

% update project
project.lastmodified = datestr(now);
if isfield(project,'idlist')
    [vals,ind]=intersect(project.idlist,idlist);
	project.idlist(ind) = [];
    if isfield(project,'arealist')
        project.arealist(ind) = [];
    end
    
	save([projectpath '\' project.name],'project');
else
    warning('can''t delete id from list');
end

return
