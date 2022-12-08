function project = ProjectCreate(projectname,desc,modelfun,parnames)
% 
% ProjectCreate - Creates a new analysis project
%
% project = ProjectCreate('projectname')
%  creates a new analysis projecti.
%  DIRS.analysis and DIRS.user must be defined.
%
% project = ProjectCreate('projectname','description')
%  attaches a description
% project = ProjectCreate('projectname','description',modelfun,parlist)
%  attaches model
%
% parlist is a cell array short parameter description e.g. {'Rmax','c50'}
%
% project = ProjectCreate('projectname','description','modelname')
%  relies on function 'modelname' to provide 'parlist'. see testmodel.m
% 

global DIRS

if nargin>2
    model = func2str(modelfun);
end

switch nargin
case 1
    desc='';model='';parnames={};
case 2
    model='';parnames={};
case 3
    parnames={};
end

projectpath = fullfile(DIRS.analysis,projectname);
projectfilename = fullfile(projectpath,projectname);

% does project already exist?
status = exist(projectpath);

if ~status% does not exist
	if ~mkdir(DIRS.analysis,projectname)
        project = [];
        str=sprintf('Cannot create project directory ''%s''.',projectpath);
        error(str)
	end
end

odir = pwd;
cd(projectpath);
status = exist([projectfilename '.mat']);
cd(odir);

if status% does exist
    project = ProjectLoad(projectname);
    str = sprintf('Project ''%s'' already exists. Loaded.',projectname);
    warning(str);
else     
    % project information
    project.name = projectname;
    project.author = DIRS.user;
    project.date =  datestr(now);
    project.lastmodified  =  datestr(now);
    project.description = desc;

    % add code to project
    mkdir(projectpath,'tools');
    % add metadata that informs the user about how the analyses files were
    % generated
    mkdir(projectpath,'metadata');
    
    project.model = '';
    project.parnames = {};
    if exist(model)==2 % model is a single m-file
        modelpath = which(model);
		disp(['model m-file: ' modelpath]);
        toolspath = [projectpath '\tools'];
    	status = copyfile(modelpath,[toolspath '\' model]);
        if status
            project.model = ['tools\' model];
        else
            disp('model file not copied successfully');
        end
        if length(parnames);
            project.parnames = parnames;
		else
            project.parnames = eval(model);
        end
    
    else % no model provided 
        disp('no valid model source code provided');
        disp('tools directory will be empty');
	end

    project.idlist = {};
    project.arealist = {};
    
    save(projectfilename,'project');
end

return

% 	if exist(model)==7     % model is a directory
%         content = what(model);
%         disp(['model directory: ' content.path]);
%         if any(strcmp(content.m,'model.m')) % model.m is included in directory
%             copydir(model,[projectpath '\' tools]);
%             project.model = 'tools/model';
%         else % model.m is not included
%             warning('directory ignored because does not include model.m');
% 		end
