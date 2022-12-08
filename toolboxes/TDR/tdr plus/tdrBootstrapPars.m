function bootpars = tdrBootstrapPars(data,bootpars)

% Check if data is sequentially or simultaneously recorded
if isfield(data,'unit') && ~isempty(data.unit)
    
    %--- Sequential recordings ---
    
    % Number of units
    nun = length(data.unit);
    
    % Number of time samples
    npt = length(data.time);
    
    % Initialize
    bootpars.iboot = cell(1,bootpars.nboot+1);
    
    % Loop over units
    for iun = 1:nun
        
        % Number of trials
        ntr = size(data.unit(iun).response,1);
        
        % Resample trials
        iboot = zeros(bootpars.nboot+1,ntr,'int32');
        iboot(1,:) = 1:ntr;
        iboot(2:end,:) = bootsamples(ntr,bootpars.nboot,[]);
        
        % Keep
        bootpars.iboot{iun} = iboot;
        
    end
    
else
    
    %--- Simultaneous recordings ---
    
    % Dimensions
    [nun,npt,ntr] = size(data.response);
    
    % Bootstrap trials
    iboot = zeros(bootpars.nboot+1,ntr);
    iboot(1,:) = 1:ntr;
    iboot(2:end,:) = bootsamples(ntr,bootpars.nboot,[]);
    
    % Keep
    bootpars.iboot = iboot;
    
end

