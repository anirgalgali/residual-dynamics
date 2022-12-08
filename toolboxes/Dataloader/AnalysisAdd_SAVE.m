function bool = AnalysisAdd(projectname, unitid, data)
% function bool = AnalysisAdd(projectname, unitid, data)
%
% AnalysisAdd  - Add new analysis (cell) to an existing project
%
% analysis = AnalysisAdd(projectname, unit, data);
%  adds a new analysis for specified unit. only one analysis per unit
%  allowed.
%

bool = 0;
global DIRS;
project = ProjectLoad(projectname);

if ~strcmp(project.author,DIRS.user)
    warning('you do not own this project. analysis create aborted.');
    return
end

projectpath = fullfile(DIRS.analysis,project.name);

analysis = [];
analysis.project = project.name;

analysisfn = [analysisfullpathname(projectname,unitid) '.mat'];

if exist(analysisfn) % analysis does not exist yet
    disp('WARNING: an analysis file for this unit id already exists');
    return;
end

analysis.date = datestr(now);
analysis.author = DIRS.user;

requiredfields = {'id','iseries','iexp','model','area'};

for ifield = 1:length(requiredfields)
    thisfield = requiredfields{ifield};
    if isfield(data,thisfield)
        analysis = setfield(analysis,thisfield,getfield(data,thisfield));
    else
        if isfield(data,'fit') & isfield(data.fit,thisfield)
            analysis = setfield(analysis,thisfield,getfield(data.fit,thisfield));
        else
            %fprintf(1,'Warning: can''t find field %s in data structure\n',thisfield);
        end
    end
end

% otherfields = setdiff(fieldnames(data),requiredfields);
% 
% for ifield = 1:length(otherfields)
%     thisfield = otherfields{ifield};
%     analysis = setfield(analysis,thisfield,getfield(data,thisfield));
% end

% Mante changed this, so that the order of the fieldnames is preserved
originalfields = fieldnames(data);

for ifield = 1:length(originalfields)
   thisfield = originalfields{ifield};
   if ~ismember(thisfield,requiredfields)
      analysis = setfield(analysis,thisfield,getfield(data,thisfield));
   end
end

AnalysisSave(projectname,unitid,analysis);

% update project
project.lastmodified = datestr(now);
if isfield(project,'idlist')
    project.idlist{end+1} = unitid;
else
    project.idlist = {unitid};
end

if isfield(project,'arealist')
if isfield(data,'area')
    project.arealist{end+1}=data.area;
else
    project.arealist{end+1}='';
end
end

save([projectpath '\' project.name],'project');

bool = 1;

return 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate percentage of explained variance 
% 
function percvar = calcpercvar(resp,predresp);

percvar = 0;

if ~isempty(resp) & ~isempty(predresp)
    if length(resp{1})==1
        vresp = cell2mat(resp);
        vpredresp = cell2mat(predresp);
        percvar = 100*(1 - mean((vresp-vpredresp).^2)/var(vresp));
    end
end
return
