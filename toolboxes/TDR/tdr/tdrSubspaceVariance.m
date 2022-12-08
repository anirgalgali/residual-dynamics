function varSUB = tdrSubspaceVariance(data_SUB,data_ORG,varpars,plotflag)
% tdrSubspaceVariance variance of responses along subspace dimensions
%
% Inputs:
%  data_SUB: population response in subspace (simultaneous, 'tdr help')
%  data_ORG: population response in full state space (simultaneous, 'tdr help')
%  varpars.time_var: index of times contributing to total variance
%  varpars.dim_sub: index of subspace dimensions to plot
%  plotflag: summary plot yes (1) or no (0)
%
% Outputs:
%  varSUB.all_nrmtot: variance(dims). 100%=total(dims,time,trial)
%  varSUB.tme_nrmtot: variance(dims,time). 100%=total(dims,time,trial)
%  varSUB.tme_nrmtme: variance(dims,time). 100%(tme)=total(dims,trial)
%
% varSUB = tdrSubspaceVariance(data_SUB,data_ORG,varpars,plotflag)

% Plot
if nargin<3
    plotflag = 0;
end

% Dimensions
[ndm npt ntr] = size(data_SUB.response);

% Parameters
if nargin<2
    varpars = [];
end
if isempty(varpars) || ~isfield(varpars,'time_var') || isempty(varpars.time_var)
    varpars.time_var = 1:npt;
end
if isempty(varpars) || ~isfield(varpars,'dim_sub') || isempty(varpars.dim_sub)
    varpars.dim_sub = 1:ndm;
end

% The response - unfolded
resp_SUB = data_SUB.response(varpars.dim_sub,:,:);
resp_ORG = data_ORG.response(:,varpars.time_var,:);

% The response - folded
Resp_SUB = resp_SUB(:,:);
Resp_ORG = resp_ORG(:,:);

% The total variance over trials and times (original state space)
totvar = nanvar(Resp_ORG,2);

% The total variance over trials (original state space)
totvar_time = nanvar(data_ORG.response,3);

% The total variance over time
refvar_time = repmat(sum(totvar_time,1),[size(resp_SUB,1) 1]);

% Output
varSUB = [];
varSUB.all_nrmtot = nanvar(Resp_SUB,2) / sum(totvar);
varSUB.tme_nrmtot = nanvar(resp_SUB,3) / sum(totvar);
varSUB.tme_nrmtme = nanvar(resp_SUB,3) ./ refvar_time;

% Plot
if plotflag
    
    figure; ha = [];
    
    % Cumulative total variance explained by each dimension
    % 100% = total variance over dimensions, times, trials
    ha(1) = subplot(2,2,1); 
    plot(100* cumsum(varSUB.all_nrmtot),'o-');
    xlabel('dimension'); ylabel('cumulative variance [%tot]');
    
    % Total variance explained by each dimension
    % 100% = total variance over dimensions, times, trials
    ha(2) = subplot(2,2,3);
    plot(log10(100* varSUB.all_nrmtot),'o-');
    xlabel('dimension'); ylabel('variance explained (log10[%tot])');
    
    % Total variance explained by each dimension - at each time
    % 100% = total variance over dimensions, times, trials
    ha(3) = subplot(2,2,2);
    hp = plot(data_ORG.time,100* varSUB.tme_nrmtot','-');
    lc = linecolors(length(hp),'jet');
    for ip = 1:length(hp)
        set(hp(ip),'color',lc(ip,:));
    end
    xlabel('time (s)'); ylabel('variance explained [%tot]');
    
    % Total variance explained by each dimension - at each time
    % 100% = total variance over dimensions and trials - at each time
    ha(4) = subplot(2,2,4); hold on;
    hp = plot(data_ORG.time,sum(100* varSUB.tme_nrmtme',2),'k--');
    hp = plot(data_ORG.time,100* varSUB.tme_nrmtme','-');
    lc = linecolors(length(hp),'jet');
    for ip = 1:length(hp)
        set(hp(ip),'color',lc(ip,:));
    end
    xlabel('time (s)'); ylabel('variance explained [%tot(time)]');
    
    set(ha,'plotbox',[1 1 1],'xlim',[-inf inf]);
    
end

