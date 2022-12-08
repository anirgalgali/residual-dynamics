function [idlist, nicelist] = ProjectGetIdList(projectname)

% Load the idlist of a project
%
% [idlist, nicelist] = ProjectGetIdList(projectname)
%
% 2003-03-17 VM made it

project = ProjectLoad(projectname);

if isfield(project,'idlist') & ~isempty(project.idlist)
   idlist = project.idlist;
   ncells = length(idlist);
   
   point = '.';
   mypoint = point + 0;
   
   for icell = 1:ncells
      myid = idlist{icell} + 0;
      ipoint = find(myid == mypoint);
      
      nicelist(icell).animal = char(myid(1:ipoint(1)-1));
      nicelist(icell).ichan = str2num(char(myid(ipoint(1)+1:ipoint(2)-1)));
      nicelist(icell).icell = str2num(char(myid(ipoint(2)+1:end)));
      nicelist(icell).id = idlist{icell};
      analysis = AnalysisLoad(projectname,idlist{icell});
      nicelist(icell).iseries = analysis.iseries;
      nicelist(icell).iexp = analysis.iexp;
   end
else
   disp(sprintf('There is no idlist in project %s',projectname));
   idlist = [];
   nicelist = [];
end

   
   