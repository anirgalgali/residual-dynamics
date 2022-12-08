function analyses = AnalysisLoad_v2(projectname,arg2,flag)
% AnalysisLoad - Load an existing analysis from a project
%
%  analyses = AnalysisLoad(projectname)
%  returns all analyses stored in project
%
%  analyses = AnalysisLoad(projectname,idlist)
%  returns analyses specified by idlist
%
%  analyses = AnalysisLoad(projectname,brainarea)
%  where brainarea = 'lgn' or 'v1'. returns all analyses with field brainarea.

global DIRS;

if nargin < 3 || isempty(flag)
    flag = 'slow';
end

if strcmp(flag,'fast');
    idlist = {arg2};
else
    project = ProjectLoad_v2(projectname);
    brainarea = '';
    if nargin < 2
        idlist = project.idlist;
    else
        if ischar(arg2)
            if isempty(intersect(arg2,{'lgn','v1'}))
                idlist = {arg2};
            else
                brainarea = arg2;
                idlist = {};
            end
        elseif iscell(arg2)
            idlist = arg2;
        end
    end

    if ~isempty(brainarea)
        if isfield(project,'arealist')
            idlist = project.idlist(find(strcmp(project.arealist,brainarea)));
        else
            disp('warning: no field arealist defined in project structure');
        end
    end
end

projectpath = fullfile(DIRS.analysis,projectname);

ianal = 1;
for index =1:length(idlist)
    thisid = idlist{index};

%     analysisfn = analysisfullpathname(projectname,thisid);
    analysisfn = fullfile(projectpath,thisid);
    if exist([analysisfn '.mat'])
        load(analysisfn);
        if size(analysis) == 1;
            analyses(ianal) = analysis;
        else
            analyses(ianal,:) = analysis;
        end
        ianal=ianal+1;
    else
        fprintf(1,'%s cannot find file %s\n',thisid,analysisfn);
    end
end

if ~exist('analyses');analyses=[];end;

return