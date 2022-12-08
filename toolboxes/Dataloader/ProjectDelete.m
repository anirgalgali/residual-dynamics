function bool = ProjectDelete(projectnames,options);
% ProjectDelete - Delete an existing analysis project. 
%  bool = ProjectDelete(projectname);
%
% You must the project owner
% options(1) = 1: batch mode, no prompts

nopt =1 ;

if nargin<2;options=zeros(1,nopt);end;

global DIRS

if ischar(projectnames)
    projectnames = {projectnames};
end

for iproj = 1:length(projectnames)
    projectname = projectnames{iproj};
    
	project = ProjectLoad(projectname);
	
	if ~isstruct(project)
        project.name = projectname;
        project.author = DIRS.user;
	end
	
	projectpath = [DIRS.analysis '\' project.name];
	
	bool = 0;
	if exist(projectpath)==7
		if strcmp(project.author,DIRS.user)
            if ~options(1)
                prompt = sprintf('Are you sure you want to delete project %s from directory %s (yes or no)? ',project.name, DIRS.analysis);
                str = input(prompt,'s');
            else
                str = 'yes';
            end
            if strcmpi(str,'yes');
                rmdir(projectpath,'s'); % requires MATLAB R14
                % deltree(projectpath,'quiet'); % chorale toolbox
                disp(sprintf('deleted project %s',project.name));
                bool = 1;
            else
                disp('delete aborted.');
            end    
		else
            warning('you do not own this project. abort.');
		end
	else
		disp('project directory does not exist');
	end
end

return