function dataTB = tdrSurrogateData(dataT,dataC,surrogateType,nboot)
% tdrSurrogateData create surrogate data
%
% Inputs:
%  dataT: trial-by-trial responses. Provides the NOISE of the responses
%  dataC: condition-averaged responses. Provides the MEAN of the responses
%  surrogateType: type of surrogate data
%  nboot: number of repetitions of the original data in the surrogate
%
% Output:
%  dataTB: trial-by-trial surrogate data. MEAN + NOISE.
%
% dataTB = tdrSurrogateData(dataT,dataC,surrogateType,nboot)

if nargin<3 || isempty(surrogateType)
    surrogateType = 'AveragePlusNoise';
end
if nargin<4 || isempty(nboot)
    nboot = 1;
end

switch surrogateType
    case 'AveragePlusNoise'
        
        % Check if data is sequentially or simultaneously recorded
        if isfield(dataT,'unit') && ~isempty(dataT.unit)
            
            %--- Sequential recordings ---
            % Number of units
            nun = length(dataT.unit);
            
            % Initialize
            dataTB = dataT;
            
            % Loop over units
            for iun = 1:nun
                
                % Trial data and condition averages
                dT = dataT.unit(iun);
                dC = dataC;
                
                % Noise response
                noise = dT.response;
                
                % Sizes
                [ntr,npt] = size(dT.response);
                
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
                response = zeros([npt ntr nboot]);
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
                    response(:,logical(jj_cnd(icd,:)),:) = ...
                        repmat(dC.response(iun,:,icd),[1 1 sum(jj_cnd(icd,:)) nboot]) + ...
                        noise_boot;
                end
                
                % Keep what you need
                for ivr = 1:nvr
                    dataTB.task_variable.(variable_name{ivr}) = task_variable.(variable_name{ivr})(:);
                    dataTB.task_index.(variable_name{ivr}) = task_index.(variable_name{ivr})(:);
                end
                dataTB.unit(iun).response = response(:,:);
                
            end
            
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


