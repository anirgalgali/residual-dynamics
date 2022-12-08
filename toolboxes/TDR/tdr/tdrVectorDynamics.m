function [autoCorr,crosCorr,h] = tdrVectorDynamics(coef,timepars,plotflag)
% tdrVectorDynamics Auto and cross correlograms of regression coefficients
%
% Inputs:
%  coef: regression coefficients
%     .name: regressor names {nrg 1}
%     .response: coefficients [ndm npt nrg]
%     .time: time axis [1 npt]
%     .dimension: dimension labels {ndm 1}
%  timepars.name: regressors to plot
%  plotflag: summary plot yes (1) or no (0)
%
% Outputs:
%  autoCorr(i): autocorrelogram for regressor i
%     .rho: correlation coefficient
%     .pval: p-value
%  crosCorr(i,j): crosscorrelograms between regressors i and j
%     .rho: correlation coefficient
%     .pval: p-value
%  h: figure, axes, and image handles
%
% [autoCorr,crosCorr,h] = tdrVectorDynamics(coef,timepars,plotflag)


if nargin<3 || isempty(plotflag)
    plotflag = 0;
end

% Regressors
nrg = length(timepars.name);

% Initialize
h.f1 = []; h.a1 = []; h.i1 = [];
h.f2 = []; h.a2 = []; h.i2 = [];

% Autocorrelation
autoCorr = [];
if plotflag
    h.f1 = figure; h.a1 = zeros(1,nrg);
    nrow = ceil(sqrt(nrg));
    ncol = ceil(nrg/nrow);
end
for irg = 1:nrg
    % Find regressor to plot
    imatch = find(strcmp(coef.name,timepars.name{irg}));
    
    % Autocorrelation
    [rho,pval] = corr(coef.response(:,:,imatch));
    
    % Keep
    autoCorr(irg).rho  = rho;
    autoCorr(irg).pval = pval;
    
    % Plot
    if plotflag
        h.a1(irg) = subplot(nrow,ncol,irg);
        h.i1(irg) = imagesc(coef.time,coef.time,abs(rho));
        caxis([0 1]);
        
        % Label
        ht=title(timepars.name{irg});
        set(ht,'interpreter','none');
        xlabel('time (s)'); ylabel('time (s)');
        set(gca,'plotbox',[1 1 1],'ydir','normal');
    end
end

% Cross correlation
crosCorr = [];
if plotflag
    h.f2 = figure;
    h.a2 = zeros(nrg,nrg);
end
for irg1 = 1:nrg
    for irg2 = 1:nrg
        
        % Find regressors to correlate
        im1 = find(strcmp(coef.name,timepars.name{irg1}));
        im2 = find(strcmp(coef.name,timepars.name{irg2}));
        
        % Cross-correlation
        [rho,pval] = corr(coef.response(:,:,im1),coef.response(:,:,im2));
        
        % Keep
        crosCorr(irg1,irg2).rho  = rho;
        crosCorr(irg1,irg2).pval = pval;
        
        % Plot
        if plotflag
            h.a2(irg1,irg2) = subplot(nrg,nrg,nrg*(irg1-1)+irg2);
            h.i2(irg1,irg2) = imagesc(coef.time,coef.time,abs(rho));
            caxis([0 1]);
            
            % Label
            hx=xlabel(coef.name{im2});
            hy=ylabel(coef.name{im1});
            set([hx hy],'interpreter','none');
        end
    end
end
if plotflag
    set(h.a2,'plotbox',[1 1 1],'ydir','normal');
    set(h.a2(1:end-1,:),'xticklabel',[]);
    set(h.a2(:,2:end),'yticklabel',[]);
end



