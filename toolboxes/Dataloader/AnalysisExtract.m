function [results] = AnalysisExtract(projname,idlist,fields,index);
% [results] = AnalysisExtract(projname,idlist,fields,index);

if length(idlist)
    array = AnalysisLoad(projname,idlist);
else
    array = AnalysisLoad(projname);
end

% if iscell(field) % legacy code
%     field = field{:};
% end

for ielement = 1:length(array)
    fieldval = getfield(array(ielement),fields{:});
    results(ielement)=fieldval(index);
end

return;