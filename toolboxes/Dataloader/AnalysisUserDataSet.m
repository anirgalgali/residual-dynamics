function bool = AnalysisUserDataSet(projectname,unit,varargin)
%
% function bool = AnalysisUserDataSet(projectname,unit,key1,val1,...)
%
%

global DIRS;
bool = 0;

if round(nargin/2) ~= nargin/2
    error('nargin should be even');
end

project = ProjectLoad(projectname);

if ~strcmp(project.author,DIRS.user)
    warning('you do not own this project. user data add aborted.');
    return
end

analysis = AnalysisLoad(projectname,unit);
projectpath = [DIRS.analysis '\' analysis.project];

iarg = 1;
while iarg < nargin - 2
    key = varargin{iarg};
    val =  varargin{iarg+1};
    if ischar(key)
        if ~isempty(analysis.userdata) % already stored data ?
            index = find(strcmp(analysis.userdata(1,:),key));
        else % no data 
            index = [];
        end
        
        if isempty(index)
            [ij,ik]=size(analysis.userdata);
            analysis.userdata{1,ik+1} = key;
            analysis.userdata{2,ik+1} = val;
            bool = 1;
        else
            warning(sprintf('key %s already exist',key));
        end
            
    else
        warning('user data key must be a string');
    end
    iarg = iarg +2;
end

analysis.lastmodified = datestr(now);
AnalysisSave(projectname,unit,analysis);

return