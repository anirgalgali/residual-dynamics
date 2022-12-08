function aname = AnalysisName(varargin)
% ANALYSISNAME make a name from arbitrary components

npts = length(varargin);
aname = [];
for ipts = 1:npts
    aname = [aname varargin{ipts}];
    if ipts < npts
        aname = [aname '_'];
    end
end

return

aname = AnalysisName('this','is','a','test');


    

