function bool = AnalysisSave(projectname,unitid,analysis)
%
% AnalysisSave
%
% bool = AnalysisSave(projectname,unit,analysis)
%

global DIRS;
projectpath = fullfile(DIRS.analysis,projectname);
pathname = fullfile(projectpath,unitid);
% analysisfn = analysisfullpathname(projectname,unitid);
save(pathname,'analysis','-v7.3');
return
