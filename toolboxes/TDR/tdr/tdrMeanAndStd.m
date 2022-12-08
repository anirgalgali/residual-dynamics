function [ravg,rstd] = tdrMeanAndStd(data,avgpars)
% tdrMeanAndStd compute mean and std of responses
%
% Inputs:
%  data: population response (sequential or simultaneous, see 'help tdr')
%  avgpars.trial: logical indeces of trials/conditions to use [ntr 1] (def: all)
%  avgpars.time: logical indeces of times to use [npt 1] (def: all)
%
% Outputs:
%  ravg: mean of responses across trials and times [nun 1]
%  rstd: STD of responses across trials and times [nun 1]
%
% [ravg,rstd] = tdrMeanAndStd(data,avgpars)

% Inputs
if nargin < 2
    avgpars = [];
end

% Time indeces
if isempty(avgpars) || ~isfield(avgpars,'time') || isempty(avgpars.time)
    timeSet = 1;
else
    timeSet = 0;
end
% Trial indeces
if isempty(avgpars) || ~isfield(avgpars,'trial') || isempty(avgpars.trial)
    trialSet = 1;
else
    trialSet = 0;
end


% Number of time samples
npt = length(data.time);
    
% Check if data is sequentially or simultaneously recorded
if isfield(data,'unit') && ~isempty(data.unit)
    
    
    %--- Sequential recordings ---    
    % Number of units
    nun = length(data.unit);
    
    % Initialize
    ravg = zeros(nun,1);
    rstd = zeros(nun,1);
    
    % Time indeces
    if timeSet
        avgpars.time = true([npt 1]);
    end
    
    % Loop over units
    for iun = 1:nun
        
        % Number of trials
        ntr = size(data.unit(iun).response,1);
        
        % Trial indeces
        if trialSet
            avgpars.trial = true([ntr 1]);
        end
        
        % Responses to use
        response = data.unit(iun).response(avgpars.trial,avgpars.time);

        % Mean across trials and times
        ravg(iun) = nanmean(response(:));
        
        % Standard deviation across trials and times
        rstd(iun) = nanstd(response(:));
        
    end
    
else
    
    
    %--- Simultaneous recordings --- 
    % Dimensions
    [nun npt ntr] = size(data.response);
    
    % Time indeces
    if timeSet
        avgpars.time = true([npt 1]);
    end
    
    % Trial indeces
    if trialSet
        avgpars.trial = true([ntr 1]);
    end
    
    % Responses to use
    response = data.response(:,avgpars.time,avgpars.trial);
    
    % Collapse time and conditions/trials
    rcollapse = response(:,:);
    
    % Mean across trials/conditions and times
    ravg = nanmean(rcollapse,2);
    
    % Mean across trials/conditions and times
    rstd = nanstd(rcollapse,2);
    
end



        
        


