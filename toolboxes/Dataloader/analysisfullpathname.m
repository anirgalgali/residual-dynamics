function pathname = analysisfullpathname(projectname,unitid)
%
% pathname = analysisfullpathname(projectname,unit)
%
% 2003 VB made it
% 2003-3-31 VM made it be happy even without id

global DIRS
projectpath = fullfile(DIRS.analysis,projectname);

unitid(find(unitid == '.'))='_';
filename = [unitid '_' projectname];
pathname = fullfile(projectpath,filename);

return;
