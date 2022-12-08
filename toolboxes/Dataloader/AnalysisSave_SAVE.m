function bool = AnalysisSave(projectname,unitid,analysis)
%
% AnalysisSave
%
% bool = AnalysisSave(projectname,unit,analysis)
%

global DIRS;
projectpath = [DIRS.analysis '\' analysis.project];
analysisfn = analysisfullpathname(projectname,unitid);
save(analysisfn,'analysis');
return
