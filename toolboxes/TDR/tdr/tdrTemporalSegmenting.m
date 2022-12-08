function data_tmp_segmented = tdrTemporalSegmenting(data,tmppars)
% Extracting temporal segments of trial-by-trial or condition-averaged responses
%
% Inputs:
%  data: population response (simultaneous, see 'help tdr')
%  tmppars.idxs: indices of times that need to be extracted
%
% Outputs:
%  data_tmp_segmented: temporally segmented data, same format as data.
%
% data_tmp_segmented = tdrTemporalSegmenting(data,tmppars)

% Initialize
data_tmp_segmented = data;
[nun,npt,ntr] = size(data.response);

if nargin < 2
    tmppars = [];
end

% Time indeces
if isempty(tmppars) || ~isfield(tmppars,'idxs') || isempty(tmppars.idxs)
    timeIdxs = true(npt,1);
else
    timeIdxs = tmppars.idxs;
end


   data_tmp_segmented.response = data.response(:,timeIdxs,:);
   data_tmp_segmented.time = data.time(timeIdxs);

end