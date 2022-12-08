function dataTB = tdrSurrogateData(dataT,dataC,pars)
% tdrSurrogateData create surrogate data
%
% Inputs:
%  dataT: trial-by-trial responses. Provides the NOISE of the responses
%  dataC: condition-averaged responses. Provides the MEAN of the responses
%  pars.surrogateType: type of surrogate data
%   - AveragePlusNoise: only works for simultaneous data. Does assume that
%      each trial is only part of one condition?
%   - AveragePlusPoisson: adds Poisson noise to the condition averages
%  pars.nboot: number of repetitions of the original data in the surrogate
%              if NaN nboot based on original trial count
%  pars.ravg: the avergage firing rate for each neuron [nun 1]. Def NaN
%  pars.rstd: the standard deviation of firing for each neuron [nun 1]. Def NaN
%  pars.Fs: sampling frequency (required for AveragePlusPoisson)
%
% Output:
%  dataTB: trial-by-trial surrogate data. MEAN + NOISE.
%
% dataTB = tdrSurrogateData(dataT,dataC,surrogateType,nboot)

if nargin<3 || isempty(pars)
    pars.surrogateType = 'AveragePlusNoise';
    pars.nboot = 1;
end    
   
if ~isfield(pars,'surrogateType') || isempty(pars.surrogateType)
    pars.surrogateType = 'AveragePlusNoise';
    disp('WARNING: only simultaneous recordings');
end
if ~isfield(pars,'nboot') || isempty(pars.nboot)
    pars.nboot = 1;
end
if ~isfield(pars,'ravg') || isempty(pars.ravg)
    pars.ravg = NaN;
end
if ~isfield(pars,'rstd') || isempty(pars.rstd)
    pars.rstd = NaN;
end
if ~isfield(pars,'ngain') || isempty(pars.ngain)
    pars.ngain = 1;
end

switch pars.surrogateType
    
    case 'AveragePlusPoisson'
        
        % Check if data is sequentially or simultaneously recorded
        if isfield(dataT,'unit') && ~isempty(dataT.unit)
            
            %--- Sequential recordings ---
            nun = length(dataT.unit);
            
            % Number of conditions
            ncd = size(dataC.response,3);
            
            % Task variables
            varname = fieldnames(dataC.task_variable);
            nvr = length(varname);
            
            % Initialize
            dataTB = dataT;
            
            % Loop over units
            for iun = 1:nun
                
                % Unit data
                dT = dataT.unit(iun);
                dC = dataC;
                
                % Dimensions
                [ntr, npt] = size(dT.response);
                
                % New number of trials
                if isnan(pars.nboot)
                    nbt = ceil(ntr/ncd);
                else
                    nbt = pars.nboot;
                end
                
                % Initialize
                dS = struct('response',NaN(nbt*ncd,npt));
                newname = cell(size(varname));
                for ivr = 1:nvr
                    % Variable name containing regression values
                    newname{ivr} = sprintf('%s_srg',varname{ivr});
                    
                    % Initialize
                    dS.task_variable.(varname{ivr}) = NaN(nbt*ncd,1);
                    dS.task_variable.(newname{ivr}) = NaN(nbt*ncd,1);
                    dS.task_index.(varname{ivr}) = NaN(nbt*ncd,1);
                    dS.task_index.(newname{ivr}) = NaN(nbt*ncd,1);
                end
                
                % Raw response
                rnrm = squeeze(dC.response(iun,:,:));
                
                % Convert to firing rates
                if isnan(pars.ravg) & isnan(pars.rstd)
                    % Mean and Std
                    ravg = nanmean(rnrm(:)-min(rnrm(:))) - min(rnrm(:));
                    rstd = 1;
                    
                    % Convert
                    rraw = rnrm*rstd + ravg;
                    
                elseif length(pars.ravg)==1 && length(pars.rstd)==1
                    % Mean and Std
                    ravg = pars.ravg;
                    rstd = pars.rstd;
                    
                    % Convert
                    rraw = rnrm*rstd + ravg;
                    
                elseif length(pars.ravg)==nun && length(pars.rstd)==nun
                    % Mean and Std
                    ravg = pars.ravg(iun);
                    rstd = pars.rstd(iun);
                    
                    % Convert
                    rraw = rnrm*rstd + ravg;
                    
                else
                    error('Bad dimensions for ravg or rstd');
                    
                end
                                
                % Loop over conditions
                for icd = 1:ncd
                    % Relevant trials
                    itr = (nbt*(icd-1)+1):nbt*icd;

                    % Keep task variables
                    for ivr = 1:nvr
                        % Variables
                        if isnan(dC.task_variable.(varname{ivr})(icd))
                            % Original variable
                            dS.task_variable.(varname{ivr})(itr,1) = NaN;
                            dS.task_index.(varname{ivr})(itr,1) = ...
                                dC.task_index.(varname{ivr})(icd);
                            
                            % New variable
                            dS.task_variable.(newname{ivr})(itr,1) = ...
                                nanmean(dT.task_variable.(varname{ivr}));
                        else
                            % Original variable
                            dS.task_variable.(varname{ivr})(itr,1) = ...
                                dC.task_variable.(varname{ivr})(icd);
                            dS.task_index.(varname{ivr})(itr,1) = ...
                                dC.task_index.(varname{ivr})(icd);
                            
                            % New variable
                            dS.task_variable.(newname{ivr})(itr,1) = ...
                                dC.task_variable.(varname{ivr})(icd);
                            dS.task_index.(newname{ivr})(itr,1) = ...
                                dC.task_index.(varname{ivr})(icd);
                        end
                    end
                                        
                    % Make surrogate responses
                    if ~all(isnan(rraw(:,icd)))
                        % Make all trials
                        rtot = repmat(rraw(:,icd),[1 nbt]);
                        
                        % Rectify
                        rtot = max(0,rtot);
                        
                        % Add noise
                        rtot = rtot + pars.ngain*randn(size(rtot)).*sqrt(rtot/pars.Fs)*pars.Fs;
                        
                        % Normalize back
                        rsur = (rtot-ravg)/rstd;
                        
                        % Keep response
                        dS.response(itr,:) = rsur';
                    end
                end
                
                % Task indeces for new variables
                for ivr = 1:nvr
                    [~,~,dS.task_index.(newname{ivr})(:,1)] = ...
                        unique(dS.task_variable.(newname{ivr}));
                end
                
                % Variable names
                variable_name = fieldnames(dS.task_variable);
                
                % Variable-index pairs
                nva = length(variable_name);
                for iva = 1:nva
                    % All values
                    all_values = dS.task_variable.(variable_name{iva});
                    
                    % Good values
                    notnan = ~isnan(all_values);
                    
                    % Unique good values
                    variable_unique = unique(all_values(notnan));
                    
                    % Variable vs. index pairs
                    dS.task_variable_index.(variable_name{iva}) = [variable_unique (1:length(variable_unique))'];
                end
                
                % Keep
                dataTB.unit(iun,1).response             = dS.response;
                dataTB.unit(iun,1).task_variable        = dS.task_variable;
                dataTB.unit(iun,1).task_index           = dS.task_index;
                dataTB.unit(iun,1).task_variable_index  = dS.task_variable_index;
            end
            
        else
            
            %--- Simultaneous recordings ---
            error('not ready for simultaneous recordings');
            
        end    
        
        
    case 'AveragePlusNoise'
        
        % Check if data is sequentially or simultaneously recorded
        if isfield(data,'unit') && ~isempty(data.unit)
            
            %--- Sequential recordings ---
            error('not ready for sequential recordings');
            
        else
            
            %--- Simultaneous recordings ---
            % Trial data and condition averages
            dT = dataT;
            dC = dataC;
            
            % Noise response
            noise = dT.response;
            
            % Sizes
            [ndm,npt,ntr] = size(dT.response);
            
            % Number of conditions
            ncd = size(dC.response,3);
            
            % All the task indeces
            variable_name = fieldnames(dC.task_index);
            nvr = length(variable_name);
            
            % Initialize
            dataTB = dT;
            dataTB.response = [];
            dataTB.task_variable = [];
            dataTB.task_index = [];
            
            % Initialize new indeces and variables
            task_variable = [];
            task_index = [];
            for ivr = 1:nvr
                task_variable.(variable_name{ivr}) = zeros(ntr,nboot);
                task_index.(variable_name{ivr})    = zeros(ntr,nboot);
            end
            
            % Loop over conditions
            response = zeros([ndm npt ntr nboot]);
            for icd = 1:ncd
                
                % Trials to average
                jj_var = zeros(nvr,ntr);
                
                % Loop over variables
                for ivr = 1:nvr
                    % Find matching trials
                    if dC.task_index.(variable_name{ivr})(icd)==0
                        % Take all trials
                        jj_var(ivr,:) = true(ntr,1);
                    else
                        % Find matching trials
                        jj_var(ivr,:) = ...
                            dT.task_index.(variable_name{ivr}) == ...
                            dC.task_index.(variable_name{ivr})(icd);
                    end
                end
                
                % Fullfill constraints on all indeces
                jj_cnd(icd,:) = prod(jj_var,1);
                
                % New indeces and variables
                for ivr = 1:nvr
                    task_variable.(variable_name{ivr})(logical(jj_cnd(icd,:)),:) = ...
                        dC.task_variable.(variable_name{ivr})(icd);
                    task_index.(variable_name{ivr})(logical(jj_cnd(icd,:)),:) = ...
                        dC.task_index.(variable_name{ivr})(icd);
                end
                
                % Number of noise samples
                nsamples = size(noise,3);
                
                % Number of trials for this condition
                ncdtr = sum(jj_cnd(icd,:));
                
                % Loop over dimensions
                ndm = size(noise,1);
                for idm = 1:ndm
                    % Bootstraped noise
                    iboots = bootsamples(nsamples,nboot,idm*100+icd);
                    % Note that the noise is being shuffled across conditions
                    % but not across time (maintains the temporal
                    % autocorrelation of the noise).
                    % In principle could also shuffle trials within conditions,
                    % however that can be problematic for conditions with small
                    % number of trials.
                    iboot  = iboots(:,1:ncdtr)';
                    noise_boot = zeros([1 npt ncdtr nboot]);
                    for jb = 1:nboot
                        noise_boot(1,:,:,jb) = noise(idm,:,iboot(:,jb));
                    end
                    
                    % The responses
                    response(idm,:,logical(jj_cnd(icd,:)),:) = ...
                        repmat(dC.response(idm,:,icd),[1 1 sum(jj_cnd(icd,:)) nboot]) + ...
                        noise_boot;
                end
            end
            
            % Keep what you need
            for ivr = 1:nvr
                dataTB.task_variable.(variable_name{ivr}) = task_variable.(variable_name{ivr})(:);
                dataTB.task_index.(variable_name{ivr}) = task_index.(variable_name{ivr})(:);
            end
            dataTB.response = response(:,:,:);
            
        end
        
    otherwise
        disp('Unknown surrogate type');
end


