function [dataPC] = tdrPca(data,pcapars,plotflag)
% tdrPca Principal components analysis of population responses
%
% Inputs:
%  data: population response (simultaneous, see 'help tdr')
%  pcapars.trial_pca: trials/conditions to compute PCs
%  pcapars.trial_prj: trials/conditions to project onto PCs
%  pcapars.time_pca: times to compute PCs
%  pcapars.plot_dimensions: dimensions to plot
%  plotflag: summary plot yes (1) or no (0)
%
% Outputs:
%  dataPC: population response in PC space
%  pc2un: projection matrix from PC space to Unit space
%  varnrm: total variance explained along each PC (normalized)
%
% [dataPC,pc2un,variances] = tdrPca(data,pcapars,plotflag)

% Also has field itime (time indeces to use)

% Time to use for PCA
if ~isfield(pcapars,'time_pca') || isempty(pcapars.time_pca)
    pcapars.time_pca = true(1,size(data.response,2));
end

% Default inputs
if nargin<3
    plotflag = 0;
end

% Dimensions
[nun npt ntr] = size(data.response);
dt = data.time(2) - data.time(1);

% Trials to use for PCA
if nargin<2 || isempty(pcapars) || ~isfield(pcapars,'trial_pca') || ...
        isempty(pcapars.trial_pca)
    % Use all trials
    trial_pca = 1:ntr;
else
    % Use specified subset of trials
    trial_pca = pcapars.trial_pca;
end

% Trials to use for projection
if nargin<2 || isempty(pcapars) || ~isfield(pcapars,'trial_prj') || ...
        isempty(pcapars.trial_prj)
    % Use all trials
    trial_prj = 1:ntr;
else
    % Use specified subset of trials
    trial_prj = pcapars.trial_prj;
end


% The responses
respUN_pca = data.response(:,pcapars.time_pca,trial_pca);
respUN_prj = data.response(:,:,trial_prj);

% Fold time and trials
RespUN_pca = respUN_pca(:,:);
RespUN_prj = respUN_prj(:,:);

% Compute PCA - select conditions
[pc2un,scores,variances] = princomp(RespUN_pca');

% Normalize variance
varnrm = variances / sum(variances);

% Project responses into PC space
% Only trials used to compute PCs (same as scores)
RespPC_pca = (pc2un'*(RespUN_pca - repmat(nanmean(RespUN_pca,2),[1 size(RespUN_pca,2)])))';
% All trials
RespPC_prj = (pc2un'*(RespUN_prj - repmat(nanmean(RespUN_pca,2),[1 size(RespUN_prj,2)])))';

% Variances
% var_pca = nanvar(RespPC_pca,1); % same as variances

% All trials in PC space - unfolded
% respPC_prj = reshape(RespPC_prj',size(respUN_prj));

% Keep responses
dataPC.response = reshape(RespPC_prj',size(respUN_prj));

% Keep task variables
varnames = fieldnames(data.task_variable);
for ivr = 1:length(varnames)
    dataPC.task_variable.(varnames{ivr}) = ...
        data.task_variable.(varnames{ivr})(trial_prj);
end

% Keep task indeces
indnames = fieldnames(data.task_index);
for idx = 1:length(indnames)
    dataPC.task_index.(indnames{idx}) = ...
        data.task_index.(indnames{idx})(trial_prj);
end

% Keep ntrial (if simultaneous recordings)
if isfield(data,'n_trial') && ~isempty(data.n_trial)
    if size(data.n_trial,1)==1
        dataPC.n_trial = data.n_trial(:,trial_prj);
    else
        % For sequential recordings n_trial is only meaningful at the level of
        % individual neurons, not PC dimensions.
        dataPC.n_trial = [];
    end
end

% Keep time
dataPC.time = data.time;
dataPC.PrinVecs = pc2un;
dataPC.PCVar = varnrm;
dataPC.responseSpeed = (dataPC.response(:,3:end,:) - dataPC.response(:,1:end-2,:))./(2*dt);
% Dimension names
dataPC.dimension = cell(nun,1);
for iun = 1:nun
    dataPC.dimension{iun,1} = sprintf('pc_%i',iun);
end

% Plot
if plotflag
    
    figure; ha = [];
    
    % Dimensions to plot
    if isempty(pcapars.plot_dimensions)
        idp = 1:20;
    else
        idp = pcapars.plot_dimensions;
    end
    
    ha(1) = subplot(1,2,1); 
    plot(100* cumsum(variances(idp) / sum(variances)),'o-');
    xlabel('principal component'); ylabel('cumulative variance [%tot]');
    
    ha(2) = subplot(1,2,2);
    plot(log10(100* variances(idp) / sum(variances)),'o-');
    xlabel('principal component'); ylabel('variance explained (log10[%tot])');
        
    set(ha,'plotbox',[1 1 1],'xlim',[-inf inf]);
end




