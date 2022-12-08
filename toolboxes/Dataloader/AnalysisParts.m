function [pcel,pstr] = AnalysisParts(aname,varargin)

% The separator
csep = double('_');

% The analysis name
anum = double(aname);

% Make sure the first and last character is a separator
if anum(1) ~= csep
    anum = [csep anum];
end
if anum(end) ~= csep
    anum = [anum csep];
end

% The separators
jsep = find(anum == csep);
nsep = length(jsep);

% Find parts
for isep = 1:nsep-1
    pcel{isep} = char(anum(jsep(isep)+1:jsep(isep+1)-1));
end

npt = length(pcel);
pstr = [];
if length(varargin) == npt
    for ipt = 1:npt
        pstr = setfield(pstr,varargin{ipt},pcel{ipt});
    end
end

return


aname = AnalysisName('this','is','a','test');

[pcel,pstr] = AnalysisParts(aname,'a1','a2','a3','a4');

