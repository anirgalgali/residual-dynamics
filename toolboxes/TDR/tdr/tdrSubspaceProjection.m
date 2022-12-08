function [data_subSUB,data_subORG] = tdrSubspaceProjection(data_fulORG,subpars,plotflag)
% tdrSubspaceProjection Project population response into subspace
%
% Inputs:
%  data_fulORG: population response in orignal space (original basis)
%  subpars.subSUB_fulORG: projection matrix from original space (original
%     basis) into subspace (subspace basis)
%  subpars.dimension: labels of subspace dimensions {ndm 1}
% 
%  Outputs:
%   data_subSUB: population response in subspace (subspace basis)
%   data_subORG: population response in subspace (original basis)
%   varnrm: total variance along subspace axes (normalized)
%
% [data_subSUB,data_subORG,varnrm] = tdrSubspaceProjection(data_fulORG,subpars,plotflag)

% Projection from original into subspace basis
if ndims(subpars.subSUB_fulORG)==2
    subSUB_fulORG = zeros([size(subpars.subSUB_fulORG,1) 1 size(subpars.subSUB_fulORG,2)]);
    subSUB_fulORG(:,1,:) = subpars.subSUB_fulORG;
elseif ndims(subpars.subSUB_fulORG)==3
    subSUB_fulORG = subpars.subSUB_fulORG;
else
    error('bad projection');
end

% Projection from subspace into original basis
subORG_subSUB = permute(subSUB_fulORG,[3 2 1]);

% Inputs
if nargin<3 || isempty(plotflag)
    plotflag = 0;
end

% Dimensions
[~,npt,ntr] = size(data_fulORG.response);

% Fixed or time-dependent projection matrix
ntt = size(subSUB_fulORG,2);
isfixed = ntt==1;

% Fold time and trials/conditions
Resp_fulORG = data_fulORG.response(:,:);

% Initialize
data_subSUB = data_fulORG;
data_subORG = data_fulORG;

% Keep ntrial (if simultaneous recordings)
if isfield(data_fulORG,'n_trial') && ~isempty(data_fulORG.n_trial)
    if size(data_fulORG.n_trial,1)==1
        data_subSUB.n_trial = data_fulORG.n_trial;
        data_subORG.n_trial = data_fulORG.n_trial;
    else
        % For sequential recordings n_trial is only meaningful at the level of
        % individual neurons, not along arbitrary state space dimensions.
        data_subSUB.n_trial = [];
        data_subORG.n_trial = [];
    end
end

% Project data into subspace and back
if isfixed
    % Respresent subspace data in subspace basis
    Resp_subSUB = squeezedim(subSUB_fulORG,2) * Resp_fulORG;
    
%     Resp_subSUB = squeeze(subSUB_fulORG,2) * Resp_fulORG;

    % Represent subspace data in original basis
    Resp_subORG = squeezedim(subORG_subSUB,2) * Resp_subSUB;
%     Resp_subORG = squeeze(subORG_subSUB,2) * Resp_subSUB;

else
    % Represent subspace data in subspace basis
    resp_subSUB = zeros([size(subSUB_fulORG,1) npt ntr]);
    for itt = 1:ntt
        resp_subSUB(:,itt,:) = squeezedim(subSUB_fulORG(:,itt,:),2) * ...
            squeezedim(data_fulORG.response(:,itt,:),2);
%           resp_subSUB(:,itt,:) = squeeze(subSUB_fulORG(:,itt,:),2) * ...
%             squeeze(data_fulORG.response(:,itt,:),2);
    end
    Resp_subSUB = resp_subSUB(:,:);
    
    % Represent subspace data in original basis
    resp_subORG = zeros(size(data_fulORG.response));
    for itt = 1:ntt
        resp_subORG(:,itt,:) = squeezedim(subORG_subSUB(:,itt,:),2) * ...
            squeezedim(resp_subSUB(:,itt,:),2);
%         resp_subORG(:,itt,:) = squeeze(subORG_subSUB(:,itt,:),2) * ...
%             squeeze(resp_subSUB(:,itt,:),2);

    end
    Resp_subORG = resp_subORG(:,:);
    
end
    
% Unfold time and trials/conditions
data_subSUB.response = reshape(Resp_subSUB,[size(Resp_subSUB,1) npt ntr]);
data_subORG.response = reshape(Resp_subORG,[size(Resp_subORG,1) npt ntr]);

% Variance along subspace axes
% varSUB = nanvar(Resp_subSUB,2)';
% 
% % Total variance
% varORG = nanvar(Resp_fulORG,2)';
% 
% % Normalize variance
% varnrm = varSUB / sum(varORG);

% The dimension names
% data_subSUB.dimension = subpars.dimension;
% data_subORG.dimension = data_fulORG.dimension;

% Plot
% if plotflag
%     
%     figure; ha = [];
%     
%     ha(1) = subplot(1,2,1); 
%     plot(100* cumsum(varSUB / sum(varORG)),'o-');
%     xlabel('subspace axis'); ylabel('cumulative variance [%tot]');
%     
%     ha(2) = subplot(1,2,2);
%     plot(log10(100* varSUB / sum(varORG)),'o-');
%     xlabel('subspace axis'); ylabel('variance explained (log10[%tot])');
%         
%     set(ha,'plotbox',[1 1 1],'xlim',[-inf inf]);
%     
% end

end
