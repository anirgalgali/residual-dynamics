function [exclude_inds,varargout] = excludeSilentUnits(data,threshold)

% This function computes the total spike count and firing rate over the
% entire duration of a trial and flags those neurons that have a mean FR 
% (computed across trials) that is lesser than that specified by threshold 
    
% INPUTS:
%
% data      - default data structure containing the response and the trial flags
%             and time event references
%               - response - N x T x K array : N - # of units, T - # of time
%                            steps, K - # of trials
% threshold - firing rate threshold that determines if a unit is classified
%            as "quiescent"
%
% OUTPUTS:
%
% exclude_inds - N x 1 logical array with '1' entrries corresponding to 
%                 quiescent units
% varargout    - cell array of variable output arguments
%                   - varargout{1} - N x K matrix of firing rate per trial         
%                     (1 time bin = T steps long)
%                   - varargout{2} - N x 1 array of mean firing rate across
%                     trials (1 time bin = T steps long)
% 
% 2018_01_18 - First written
%
% Author: Aniruddh Galgali
%
totalSpikeCountperTrial = squeeze(sum(data.response,2));
spikeFRperTrial = totalSpikeCountperTrial./(data.time(end) - data.time(1));


meanSpikeRateAcrossTrials = mean(spikeFRperTrial,2);
exclude_inds = (meanSpikeRateAcrossTrials <= threshold);

varargout{1} = spikeFRperTrial;
varargout{2}  = meanSpikeRateAcrossTrials;
    
%     data_out.response(exclude_inds,:,:) = [];
%     data_out.dimension(exclude_inds,:,:) = [];
    
end