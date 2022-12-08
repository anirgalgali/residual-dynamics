function bool = AnalysisUserDataDelete(projectname,unit,key)
%
% AnalysisUserDataDelete
%
% bool = AnalysisUserDataDelete(projectname,unit,key)
%

global DIRS;
bool = 0;

[analysis,project] =  AnalysisLoad(projectname,unit);

if ~strcmp(project.author,DIRS.user)
    warning('you do not own this project. user data add aborted.');
    return
end

index = find(strcmp(analysis.userdata(1,:),key));
if length(index)
	analysis.userdata(:,index) = [];
	AnalysisSave(projectname,unit,analysis);
    bool = 1;
end

return
