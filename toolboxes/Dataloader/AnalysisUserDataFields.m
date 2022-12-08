function fieldsarray = AnalysisUserDataFields(projectname,unit)
%
% val = AnalysisUserDataFields(projectname,unit)
%

analysis = AnalysisLoad(projectname,unit);

if ~isempty(analysis.userdata)
    fieldsarray = analysis.userdata(1,:);
else
    fieldsarray = {};
end

return
