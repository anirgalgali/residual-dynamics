function [out] = extract_session_segments(sessionResponse,sessionSegmentationParams)

% Determines the trial indices within a session at which the largest
% changes in firing rate occurs due to artefacts (e.g non-stationarities) 
% induced in the recording procedure (for example, by changing the threshold
% for spike detection on the array)
%
% INPUTS:
%
% sessionResponse           - N x K matrix of temporall averaged firing rates - 
%                            N = # of units, K = # of trials in a session.
% sessionSegmentationParams - structure containing parameters that determine segmentation
%            
%                                nChangePts -- number of change points to find for each unit  
%                                              in the first pass                 
%                                changeCountThreshold -- threshold to determine which trial
%                                                        indices have sufficient change points                                 
%                                numTrialAvgWin -- number of trials on either side of the  
%                                                  putative change point to determine relative
%                                                  change in firing rate 
%                                relChangeThreshold -- threshold for relative change in firing                                       
%                                                      rate that determines segment coundaries 
%                       
% OUTPUTS:
%
% out                        - output structure that contains the segment boundaries and 
%                              associated metrics that help determine it 
%                                
%                                 idxCptFinal -- Trial index determing the segment boundary
%                                 spikeFRRelChange -- relative change in  firing rate
%                                 around segment boundary
%
%   2018_01_18 - First written 
%   2022-04-03 - Modified
%   Author: Aniruddh Galgali

%%
for ii = 1: size(sessionResponse,1)
    [cptUnit{ii}] = findchangepts(sessionResponse(ii,:),'MaxNumChanges',sessionSegmentationParams.nChangePts);
end
out.cptUnit = cptUnit;
%% Finding counts for change points in each trial bin

[cptCounts,~] = histcounts(cell2mat(cptUnit),1:size(sessionResponse,2));
cptRetained = find(cptCounts > sessionSegmentationParams.changeCountThreshold);

out.cptCounts = cptCounts;
out.cptFirstPass = cptRetained;
%% if change points are too close by (within a 10 trial wondow) - choose best

if(~isempty(cptRetained))
    
    for jj = 1: length(cptRetained)
        
        idxsCptWin = cptRetained(jj)-10 : cptRetained(jj)+10;
        if(idxsCptWin(end) > length(cptCounts))
            idxsCptWin = idxsCptWin(idxsCptWin <= length(cptCounts));
        elseif(idxsCptWin(1) < 1)
            idxsCptWin = idxsCptWin(idxsCptWin >= 1);
        end
        [~,idx_max] = max(cptCounts(idxsCptWin));
        idxCptPutative(jj) = idxsCptWin(idx_max);
        
    end
    
    idxCptPutative = unique(idxCptPutative);
else
    idxCptPutative = [];
end
out.idxCptPutative = idxCptPutative;

%% refine change points by evaluating relative change in mean in a pre/post window

if(~isempty(idxCptPutative))
    
    numTrialAvgWin = sessionSegmentationParams.numTrialAvgWin;
    
    for ii = 1:length(idxCptPutative)
        
        idxsPre = idxCptPutative(ii)-1:-1:idxCptPutative(ii)-numTrialAvgWin;
        idxsPost = idxCptPutative(ii)+1:idxCptPutative(ii)+numTrialAvgWin;
        
        if(idxsPre(end) < 1)
            idxsPre = idxsPre(idxsPre >=1);
        end
        
        if(idxsPost(end) > size(sessionResponse,2))
            idxsPost = idxsPost(idxsPost <= size(sessionResponse,2));
        end
        
        if(length(idxsPre) < length(idxsPost))
            idxsPost = idxsPost(1:length(idxsPre));
        elseif(length(idxsPre) > length(idxsPost))
            idxsPre = idxsPre(length(idxsPre) - length(idxsPost)+1:end);
        end
        
        for jj = 1:size(sessionResponse,1)
            
            spikeFRFwd = sessionResponse(jj,idxsPost);
            spikeFRBkw = sessionResponse(jj,idxsPre);
            spikeFRRelChange(jj,ii) = abs((mean(spikeFRFwd,2) - mean(spikeFRBkw,2))./mean(spikeFRBkw,2));
            
        end
    end
    
    spikeFRRelChangeMedian = nanmedian(spikeFRRelChange,1);
    idxCptFinal = idxCptPutative(spikeFRRelChangeMedian > sessionSegmentationParams.thresholdRelMeanChange);
else
    idxCptFinal =[];
    spikeFRRelChange = [];
end
out.idxCptFinal = idxCptFinal;
out.spikeFRRelChange = spikeFRRelChange;


idxSegBreaks =  [1 out.idxCptFinal size(sessionResponse,2)];
idxsAllTrials = 1:size(sessionResponse,2);
for ii = 1:length(idxSegBreaks)- 1
     
     if(ii == length(idxSegBreaks) - 1)
         out.idxsSegment{ii} = idxsAllTrials >= idxSegBreaks(ii) & idxsAllTrials <= idxSegBreaks(ii+1);
         out.trialNumSegment{ii} = [idxSegBreaks(ii) idxSegBreaks(ii+1)];
     else
         out.idxsSegment{ii} = idxsAllTrials >= idxSegBreaks(ii) & idxsAllTrials < idxSegBreaks(ii+1);  
         out.trialNumSegment{ii} = [idxSegBreaks(ii) idxSegBreaks(ii+1)-1];
     end
       
 end

end