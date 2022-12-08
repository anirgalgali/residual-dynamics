% analysis base
global DIRS
DIRS.analysis = 'z:\analysis';
DIRS.user = 'vincent';




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
function val = AnalysisUserDataDelete(projectname,unit,key)
global DIRS;

project = ProjectLoad(projectname);

if ~strcmp(project.author,DIRS.user)
    warning('you do not own this project. user data add aborted.');
    return
end

analysis =  AnalysisLoad(project,unit);
index = find(strcmp(analysis.userdata(1,:),key));
analysis.userdata(:,index) = [];
AnalysisSave(projectname,unit,analysis);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pathname = analysisfullpathname(projectname,unit)
global DIRS
projectpath = [DIRS.analysis '\' projectname];
str = unit.id; str(find(str == '.'))='_';
filename = [str '_' project.name];
pathname = [projectpath '\' filename];
return;
