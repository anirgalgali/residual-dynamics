function val = AnalysisUserDataGet(projectname,unit,key)
%
% val = AnalysisUserDataGet(projectname,unit,key)
%

analysis =  AnalysisLoad(projectname,unit);


index = find(strcmp(analysis.userdata(1,:),key));
val = analysis.userdata{2,index};
        
return
