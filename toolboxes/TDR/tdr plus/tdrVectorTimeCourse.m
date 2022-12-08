function [coef_nrm,coef_spd] = tdrVectorTimeCourse(coef_time,plotflag)
% tdrVectorTimeAverage Norm and speed of time dependent regression vectors
%
% Inputs:
%  coef_time: regression coefficients 
%     .name: regressor names [nvc 1]
%     .response: coefficients [nun npt nvc]
%     .time: time axis [1 npt]
%  plotflag: summary plot yes (1) or no (0)
%
% Outputs:
%  coef_nrm: vector norm [npt nvc]
%  coef_spd: vector speed of change [npt nvc]
%
% [coef_nrm,coef_spd] = tdrVectorTimeCourse(coef_time,plotflag)


if nargin < 2 || isempty(plotflag)
    plotflag = 0;
end

% The vectors
vecname = coef_time.name;
nvc = length(vecname);
npt = length(coef_time.time);

% The norm of the vectors
coef_nrm = zeros(npt,nvc);
for ivc = 1:nvc
    for ipt = 1:npt
        coef_nrm(ipt,ivc) = norm(squeeze(coef_time.response(:,ipt,ivc)));
    end
end

% Initialize
coef_spd = zeros(npt,nvc);

% The time difference
time_dif = diff(coef_time.time);
    
% Loop over vectors
nvc = length(vecname);
for ivc = 1:nvc
    
    % The response difference
    resp_dif = diff(coef_time.response(:,:,ivc),[],2);
    
    % The response speed
    coef_spd(2:end,ivc) = (sqrt(sum(resp_dif.^2,1)) ./ time_dif)';
    
end


% PLOT a summary with the norm of the vectors and the times when they were
% averaged
if plotflag
    % Line colors
    lc = linecolors(nvc,'jet');
    
    figure;
    
    % Coefficient norm
    ha(1)=subplot(1,2,1); hold on;
    for ivc = 1:nvc

        % The raw vectors
        hp1=plot(coef_time.time,coef_nrm(:,ivc),'-');
        set(hp1,'color',lc(ivc,:));
        
        % Axis labels
        xlabel('time (s)'); ylabel('coefficient norm');
        
        % The vector name
        ht=text(coef_time.time(end),coef_nrm(end,ivc),coef_time.name{ivc});
        set(ht,'horizontalalignment','left','verticalalignment','top','interpreter','none',...
            'color',lc(ivc,:));
    end
    
    % Coefficient speed
    ha(2)=subplot(1,2,2); hold on;
    for ivc = 1:nvc
        
        % The raw vectors
        hp1=plot(coef_time.time,coef_spd(:,ivc),'-');
        set(hp1,'color',lc(ivc,:));
        
        % Axis labels
        xlabel('time (s)'); ylabel('coefficient speed');
        
        % The vector name
        ht=text(coef_time.time(end),coef_spd(end,ivc),coef_time.name{ivc});
        set(ht,'horizontalalignment','left','verticalalignment','top','interpreter','none',...
            'color',lc(ivc,:));
    end
    
end



