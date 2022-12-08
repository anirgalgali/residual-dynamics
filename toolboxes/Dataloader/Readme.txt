% analysis base
global DIRS
DIRS.analysis = 'z:\analysis';
DIRS.user = 'vincent';

function project = ProjectCreate(projectname,projectdesc)
global DIRS
projectpath = [DIRS.analysis '\' projectname];
e = exist(projectpath);
if ~e
    mkdir(projectpath);
    currentdir = pwd;
    cd(projectpath);
    mkdir tools;
    project.name = projectname;
    project.author = DIRS.user;
    project.createdon =  datestr(now);
    project.lastmodified  =  datestr(now);
    project.description = projectdesc;
    save(projectname,'project');
    cd(currentdir);
else
    project = [];
    error('Project name already exist');
end

return

%%%%

function bool = ProjectDelete(projectname);

function anal = AnalysisCreate(unit,model,boundaries);
anal.date = datestr(now);
anal.author = DIRS.user;
anal.animal = unit.animal;
anal.iseries = unit.iseries;
anal.iexp = unit.iexp;
anal.ichan = unit.ichan;

if exist(model)==2
    anal.model = model;
    parnames = 
else
    anal.model = [];
end

anal.parnames =
anal.boundaries = [];

return 

function bool = AnalysisAdd(projectname,)
	modelpath = which(model);
	if length(modelpath)
    	copyfile(modelpath,destpath);
    end

global DIRS

function