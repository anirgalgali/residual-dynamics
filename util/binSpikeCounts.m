function [data_out] = binSpikeCounts(data, binWidth, Fs, useSqrt, useRates)
% Function to bin neural spike counts using a specified bin-width
% Input 
% - data (struct) - contains the neural data and the associate task
%  variable for a single recording session.
% -binWidth(scalar) - specifies the bin width used for binning spike counts
%  (in milliseconds)
% -Fs (sacalar) - sampling frequency (inverse of step-sie) of the input data
% - useSqrt (boolean) - boolean to indicate whether to return square-root
%   transformed spike counts.
% - useRates (boolean) - to indicate whether to return firing rates instead
%   of spike counts.
% Outputs
% - data_out (struct) - output
%
% Author : Aniruddh Galgali. This function is heavily adapated from the
% GPFA toolbox developed by Byron Yu/John Cunningham.
  
yDim = size(data.response,1);

nBins = round(binWidth/((1/Fs)*1000));

T    = floor(size(data.response,2) / nBins); % number of time-points
K = size(data.response,3); % number of trials

data_out = data;
data_out.response = nan(yDim,T,K);
data_out.time = nan(1,T);

if(isfield(data_out,'time_rel'))
    data_out.time_rel = nan(1,T); 
end

if(isfield(data_out,'time_rel'))
    data_out.time_iev = nan(1,T); 
end

for kk = 1:K
    for tt = 1:T
        
        iStart = nBins * (tt-1) + 1;
        iEnd   = nBins * tt;
        data_out.response(:,tt,kk)= sum(data.response(:, iStart:iEnd,kk), 2);
%         data_out.response(:,tt,kk)= mean(data.response(:, iStart:iEnd,kk), 2);
    end
    if useSqrt
        data_out.response(:,:,kk) = sqrt(data_out.response(:,:,kk));
    end
    
    
    if(useRates)
        data_out.response(:,:,kk) = data_out.response(:,:,kk).*(1000./binWidth);
    end
    
end
for tt = 1:T
    
    iStart = nBins * (tt-1) + 1;
    iEnd   = nBins * tt;
    data_out.time(tt) = median(data.time(iStart:iEnd));
    if(isfield(data_out,'time_rel'))
        data_out.time_rel(tt) = median(data.time_rel(iStart:iEnd));
    end
    
    if(isfield(data_out,'time_iev'))
        data_out.time_iev(tt) = median(data.time_iev(iStart:iEnd));
    end
    
end
end

