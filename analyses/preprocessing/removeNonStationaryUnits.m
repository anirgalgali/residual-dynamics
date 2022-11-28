function [out,idxsNonStat] = removeNonStationaryUnits(response,nonStatParams)
% This function identifies non-stationary units in a dataset based on 
% heuristics specified in nonStatParams.
    
% INPUTS:
%
% response  - N x T x K array : N - # of units, T - # of time
%             steps, K - # of trials
% nonStatParams (struct) - contains heuristics computed based on the
%                        windowed firing rates that enable indetifying
%                        unith with non-stationarities.
% OUTPUTS:
% out (struct) - contains additional information about the firing rate
%                properties of the non-stationary units
% idxsNonStat - N x 1 logical array with '1' entrries corresponding to 
%                 non-staitonary units
% varargout    - cell array of variable output arguments
%                   - varargout{1} - N x K matrix of firing rate per trial         
%                     (1 time bin = T steps long)
%                   - varargout{2} - N x 1 array of mean firing rate across
%                     trials (1 time bin = T steps long)
% 
% 2018_01_18 - First written
%s
% Author: Aniruddh Galgali
nTrialsinSeg = size(response,2);

if(nTrialsinSeg < nonStatParams.minTrialsinSeg)
    
    warning('NOT ENOUGH TRIALS!, SKIPPING SEGMENT')
    
    out = [];
    idxsNonStat = [];
    return;
    
else
    
    
    movMean = [];
    
    idxStart = 1:nTrialsinSeg;
    idxEnd = idxStart + (nonStatParams.winLength - 1);
    idxEnd(idxEnd > nTrialsinSeg) = NaN;
    
    idxStart = downsample(idxStart,nonStatParams.winStep);
    idxEnd = downsample(idxEnd,nonStatParams.winStep);
    
    winsToKeep = nTrialsinSeg - idxStart >= nonStatParams.minTrialsToAvg;
    idxStart = idxStart(winsToKeep);
    idxEnd = idxEnd(winsToKeep);
    
    assert(length(idxStart) == length(idxEnd) , 'start and end indices do not match in length');
    
    
    for iWin = 1: length(idxStart)
        
        if(~isnan(idxEnd(iWin)))
            idxWin {iWin}  = idxStart(iWin):idxEnd(iWin);
        else
            idxWin {iWin}  = idxStart(iWin):nTrialsinSeg;
        end
        %             mov_var(:,iWin) = var(segmentResponse(:,idxWin),0,2);
        movMean(:,iWin) = mean(response(:,idxWin{iWin}),2);
        %             mov_median(:,iWin) = median(segmentResponse(:,idxWin),2);
        
    end
    
    
    movMeanDiff = (movMean(:,3:end) - movMean(:,1:end-2))/2;
    
    movMeanDiff_mean = mean(movMeanDiff(:));
    movMeanDiff_std = std(movMeanDiff(:));
    
    thresholdMask_mean = movMeanDiff >= movMeanDiff_mean + 3*movMeanDiff_std ...
        | movMeanDiff <= movMeanDiff_mean - 3*movMeanDiff_std ;
    
    thresholdMask_percentile = bsxfun(@ge,abs(movMeanDiff),prctile(abs(movMeanDiff),...
        nonStatParams.percentileThreshold,2));
    
    
    switch nonStatParams.thresholdToUse
        
        case 'mean_std'
            
            threshCrossings = thresholdMask_mean;
            
        case 'percentile'
            
            threshCrossings = thresholdMask_percentile;
            
            
    end
    
    responseRelChange = zeros(size(threshCrossings));
    
    for cc = 1: size(threshCrossings,1)
        
        idxWinCrossings = idxWin(threshCrossings(cc,:) == 1);
        idxWinCrossings_locs =  find(threshCrossings(cc,:) == 1);
        kk = 1;
        
        for jj = 1:length(idxWinCrossings)
            
            idxWinCentre = round(median(idxWinCrossings{jj}));
            
            if(idxWinCentre <= round(nonStatParams.winLength/2))
                continue;
            elseif(idxWinCentre >= size(response,2) - round(nonStatParams.winLength/2))
                break;
            else
                
                idxBkw = idxWinCentre - 1 : -1 : idxWinCentre - nonStatParams.nTrialsLookAhead;
                idxFwd = idxWinCentre + 1: idxWinCentre + nonStatParams.nTrialsLookAhead;
                
                if(idxBkw(end) < 1)
                    idxBkw = idxBkw(idxBkw >=1);
                end
                
                if(idxFwd(end) > size(response,2))
                    idxFwd = idxFwd(idxFwd <= size(response,2));
                end
                
                if(length(idxBkw) < length(idxFwd))
                    idxFwd = idxFwd(1:length(idxBkw));
                elseif(length(idxBkw) > length(idxFwd))
                    idxBkw = idxBkw(length(idxBkw) - length(idxFwd)+1:end);
                end
                
                
                responseBkw = mean(response(cc,idxBkw));
                responseFwd = mean(response(cc,idxFwd));
                
                responseRelChange(cc,idxWinCrossings_locs(jj)) = abs((responseFwd - responseBkw)/responseBkw);
                
                kk = kk +1;
                
            end
            
        end
        
    end
    
    
    idxsNonStat = any(responseRelChange >= nonStatParams.relChangeThreshold,2);
    
    
    
end


out.responseRelChange = responseRelChange;
out.thresholdCrossings = threshCrossings;
out.movMean = movMean;
out.movMeanDiff = movMeanDiff;
out.movMeanDiff_mean = movMeanDiff_mean;
out.movMeanDiff_std = movMeanDiff_std;



end

