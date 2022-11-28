function [idxLowChoiceSel,varargout] = computeChoiceSelectivity(data,threshold,numTimeSegments)

%{ Computes a choice selectivity index for each unit
% Input
% - data(struct) - containing responses and associated task/trial variables
% - threshold(scalar) - threshold that determines choice selectivity
% - numTimeSegments (scalar) - number of time-bins to split the time-course
%   into while computing selectivity.
% Output
% - idxLowChoiceSel (array,logical) - of size n_units x 1, 1s indicate
%   units with low choice selectivity.
% - varargout
% Author : Aniruddh Galgali (Jan 2018)
%}

% Rebins the coherencies for each trial into a 'low, 'medium'
% and 'high' categories
dt = data.time(2) - data.time(1);
% Splits the entire time course into a specified number of time segments
segTimes = linspace(data.time(1),data.time(end),numTimeSegments + 1);
idxLowChoiceSel = true(size(data.response,1),1);

idxChoice1Trials = data.task_index.targ_dir == 1;
idxChoice2Trials = data.task_index.targ_dir == 2;
% dprime = zeros(numTimeSegments,1);
% windowedChoiceSelectivity = zeros(numTimeSegments,1);
for ii = 1: numTimeSegments
    
    idx_t =  data.time < segTimes(ii + 1) & data.time >= segTimes(ii);
    segLength = segTimes(ii+1) - segTimes(ii);
    
    segSpikeFRChoice1 = (squeeze(sum(data.response(:,idx_t,idxChoice1Trials),2)))./segLength;
    segSpikeFRChoice2 = (squeeze(sum(data.response(:,idx_t,idxChoice2Trials),2)))./segLength;
    
    segLowChoiceIdx = abs(mean(segSpikeFRChoice1,2) - mean(segSpikeFRChoice2,2)) < threshold ;
    
    idxLowChoiceSel = idxLowChoiceSel & segLowChoiceIdx;
    
    dprime(:,ii) = (mean(segSpikeFRChoice1,2) - mean(segSpikeFRChoice2,2))./...
        (sqrt(0.5*(var(segSpikeFRChoice1,0,2) + var(segSpikeFRChoice2,0,2))));
    
    windowedChoiceSelectivity(ii).tLims = [segTimes(ii) segTimes(ii + 1)-dt];
    windowedChoiceSelectivity(ii).dprime = dprime;

end

varargout{1} = mean(dprime,2);
varargout{2} = windowedChoiceSelectivity;

end