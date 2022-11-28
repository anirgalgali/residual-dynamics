function [delayLengthsbinned,idxDelayLengths] = binDelayLengths(data,binBoundaries)

%{ Bins the delay lengths associated with trials into discrete bins
% Input
% - data(struct) - containing responses and associated task/trial variables
% - binBoundaries(array) - of size n_bins x 1 specifying the edge of the bins
% Output
% - delayLengthsbinned (array) - of size n_trials x 1, "binned"/discretized 
%   delay lengths based on the bins deined by binBoundaries. 
%
% Author : Aniruddh Galgali (May 2017)
%}

% Rebins the coherencies for each trial into a 'low, 'medium'
% and 'high' categories


delayPeriodLength = abs(data.task_event.fpoffcd - data.task_event.stimoffcd);
[~,~,idxDelayLengths] = histcounts(delayPeriodLength,binBoundaries);
uniqueBins = unique(idxDelayLengths);
delayLengthsVal = NaN(length(uniqueBins),1);

for jj = 1 : length(binBoundaries) - 1
    delayLengthsVal(jj) = 0.5*(binBoundaries(jj) + binBoundaries(jj+1));
end
delayLengthsbinned = NaN(length(delayPeriodLength),1);
% idxsdelayBinned = NaN(length(delayPeriodLength),1);
for jj = 1: length(uniqueBins)
    delayLengthsbinned(idxDelayLengths == uniqueBins(jj)) = delayLengthsVal(uniqueBins(jj));
       
end

end