function [difdata,hf,difplot] = tdrGradientField(data,difpars,plotflag)
% Compute distribution of response gradients in 2D subspace
%
% Inputs:
%  data: data structure
%  difpars:
%   .gradient_task_index: indeces to compute gradient. Each condition is
%      computed and plotted separately.
%   .average_task_index: indeces to compute condition averages
%   .average_subt: substract the condition averages
%   .average_plot: plot the condition averages
%   .dif_dim: dimensions to use [1 2]
%   .dif_dist: step size for the gradient computation. Def=1
%   .time_win: time window for grandient computation
%   .dif_type: type of gradient {'forward','backward','central'}. Def {'central'}
%   .bin_edge_x: vector of bin edges along x-direction
%   .bin_edge_y: vector of bin edges along y-direction
%   .bin_edge_n: number of bins. Only used if no bin_edge_x or bin_edge_y 
%   .bin_edge_w: range covered by bins, as fraction of response STD
%   .bin_nmin: minimum number of bin counts to compute gradient
%
% Outputs:
%  difdata.hist_response: 2D-histogram of responses
%  difdata.hist_gradient: 2D-histogram of gradients
%  difdata.hist_gradientMask: 2D-histogram of masks on the gradients
%  hf: figure handles
%  difplot: plotting parameters
%    .ax: x-grid
%    .ay: y-grid
%
% [difdata,hf] = tdrGradientField(data,difpars,plotflag)

% Check inputs
if ~isfield(difpars,'average_subt') || isempty(difpars.average_subt)
    average_subt = 0;
else
    average_subt = difpars.average_subt;
end
if ~isfield(difpars,'average_plot') || isempty(difpars.average_plot)
    average_plot = 0;
else
    average_plot = difpars.average_plot;
end
if length(difpars.dif_dim)~=2
    error('Wrong number of dimensions');
end
if isfield(difpars,'dif_dist') && ~isempty(difpars.dif_dist)
    ds = difpars.dif_dist;
else
    ds = 1;
end
if isfield(difpars,'dif_type') && ~isempty(difpars.dif_type)
    dif_type = difpars.dif_type;
else
    dif_type = 'central';
end
if isfield(difpars,'time_win') && ~isempty(difpars.time_win)
    [~,it1] = min(abs(difpars.time_win(1)-data.time));
    [~,it2] = min(abs(difpars.time_win(2)-data.time));
    itime = it1:it2;
else
    itime = 1:length(data.time);
end
if isfield(difpars,'bin_edge_n') && ~isempty(difpars.bin_edge_n)
    bin_edge_n = difpars.bin_edge_n;
else
    bin_edge_n = 20;
end
if isfield(difpars,'bin_edge_w') && ~isempty(difpars.bin_edge_w)
    bin_edge_w = difpars.bin_edge_w;
else
    bin_edge_w = 0.3;
end
if isfield(difpars,'bin_nmin') && ~isempty(difpars.bin_nmin)
    bin_nmin = difpars.bin_nmin;
else
    bin_nmin = 50;
end

% Sample duration
dt = diff(data.time(1:2));

% Number of time points
npt = length(data.time);

% Determine responses for gradient computation and condition averages
if isfield(difpars,'average_task_index') && ~isempty(difpars.average_task_index)
    % Remember
    got_average_task_index = 1;
    
    % Number of average conditions
    average_index_name = fieldnames(difpars.average_task_index);
    ncd_a = length(difpars.average_task_index.(average_index_name{1}));
    
    % Condition averages and mean subtraction
    [dataC,dataS] = tdrAverageCondition(data ,difpars.average_task_index);
    [dataCS,~]    = tdrAverageCondition(dataS,difpars.average_task_index);
    
    % Response for gradient computation
    if average_subt
        response  = dataS.response;
        responseC = dataCS.response;
    else
        response  = data.response;
        responseC = dataC.response;
    end
    
else
    % Remember
    got_average_task_index = 0;

    % Response for gradient computation
    response  = data.response;
    
    % Condition average response
    responseC = [];
end

% Conditions to average
ntr = size(response,3);
if ~isfield(difpars,'gradient_task_index') || isempty(difpars.gradient_task_index)
    % Number of conditions
    ncd_g = 1;
    
    % Trials to use
    jtrial = ones(1,ntr);
    
    % Match gradient and averaging conditions
    if got_average_task_index
        jconds = true(1,ncd_a);
    end
else
    % Number of conditions
    gradient_index_name = fieldnames(difpars.gradient_task_index);
    nindex = length(gradient_index_name);
    ncd_g = length(difpars.gradient_task_index.(gradient_index_name{1}));
    
    % Loop over conditions
    jtrial = zeros(ncd_g,ntr);
    for icd = 1:ncd_g
        
        % Loop over index names
        jtrial_index = zeros(nindex,ntr);
        for jindex = 1:nindex
            jtrial_index(jindex,:) = data.task_index.(gradient_index_name{jindex}) == ...
                difpars.gradient_task_index.(gradient_index_name{jindex})(icd);
        end
        
        % Combine indeces
        jtrial(icd,:) = prod(jtrial_index,1);
    end
    
    % Match gradient and averaging conditions
    if got_average_task_index
        jconds = false(ncd_g,ncd_a);
        for icd = 1:ncd_g
            
            jconds_index = zeros(nindex,ncd_a);
            for jindex = 1:nindex
                jconds_index(jindex,:) = ...
                    difpars.average_task_index.(gradient_index_name{jindex}) == ...
                    difpars.gradient_task_index.(gradient_index_name{jindex})(icd);
            end
            
            % Combine indeces
            jconds(icd,:) = logical(prod(jconds_index,1));
            
        end
    end
end
jtrial = logical(jtrial);

% Compute gradients
difference = NaN(size(response));
switch dif_type
    case 'central'
        % Central difference
        difference(:,ds+1:npt-ds,:) = (response(:,2*ds+1:npt,:) - response(:,1:npt-2*ds,:))/(2*dt);
    case 'forward'
        % Where the activity will be next
        difference(:,1:npt-ds,:) = (response(:,ds+1:npt,:) - response(:,1:npt-ds,:))/dt;
    case 'backward'
        % Where the activity was before
        difference(:,ds+1:npt,:) = (response(:,ds+1:npt,:) - response(:,1:npt-ds,:))/dt;
end

% Define bin edges
if isfield(difpars,'bin_edge_x') && ~isempty(difpars.bin_edge_x) && ...
        isfield(difpars,'bin_edge_y') && ~isempty(difpars.bin_edge_y)
    % Bin edges
    hx = difpars.bin_edge_x;
    hy = difpars.bin_edge_y;
else
    % Responses to consider
    rdim = response(difpars.dif_dim,itime,any(jtrial,1));
    
    % Overall STD
    rstd = std(rdim(:));
    
    % Overall mean
    % ravg = mean(rdim(:)); % VM 2015-07-02
    xavg = mean(rdim(1,:)); % VM 2015-07-02
    yavg = mean(rdim(2,:)); % VM 2015-07-02
    
    % Bin edges
    hx = xavg + linspace(-rstd*bin_edge_w,+rstd*bin_edge_w,bin_edge_n);
    hy = yavg + linspace(-rstd*bin_edge_w,+rstd*bin_edge_w,bin_edge_n);
end

% Initialize
hist_response = zeros(length(hx)-1,length(hy)-1,ncd_g);
hist_gradient = zeros(length(hx)-1,length(hy)-1,2,ncd_g);
hist_gradientMask = zeros(length(hx)-1,length(hy)-1,2,ncd_g);

% Compute vector field
for icd = 1:ncd_g
    % The trials for this condition
    rsp = response(difpars.dif_dim,itime,jtrial(icd,:));
    dif = difference(difpars.dif_dim,itime,jtrial(icd,:));
    
    % The position and gradient histograms
    [hispos,hisvfd,ax,ay] = hist2DVectorField(rsp(:,:)',dif(:,:)',hx,hy);
    
    % Do not use bins with too little samples
    [irowfew,icolfew] = find(hispos<bin_nmin);
    hismask = ones(size(hisvfd));
    for ifew = 1:length(irowfew)
        hismask(irowfew(ifew),icolfew(ifew),:) = NaN;
    end
    
    % Keep what you need
    hist_response(:,:,icd) = hispos;
    hist_gradient(:,:,:,icd) = hisvfd;
    hist_gradientMask(:,:,:,icd) = hismask;
end

% Keep plot info
difplot = [];
difplot.ax = ax;
difplot.ay = ay;

% Plot
if plotflag
    hf = zeros(1,ncd_g);
    for icd = 1:ncd_g
        % The histograms
        hpos = hist_response(:,:,icd);
        hvfd = hist_gradient(:,:,:,icd);
        hmsk = hist_gradientMask(:,:,:,icd);
        
        hf(icd)=figure; ha=[]; hi=[];
        
        ha(1)=subplot(2,2,1); hi(1)=imagesc(ax,ay,hvfd(:,:,1)'.*hmsk(:,:,1)');
        ylabel('dim2');
        hc=colorbar; axes(hc); ylabel('dim1 difference');
        
        ha(2)=subplot(2,2,2); hi(2)=imagesc(ax,ay,hvfd(:,:,2)'.*hmsk(:,:,2)');
        hc=colorbar; axes(hc); ylabel('dim2 difference');
        
        ha(3)=subplot(2,2,3); hi(3)=imagesc(ax,ay,log10(hpos(:,:)')); hold on; 
        if got_average_task_index && average_plot
            rr1 = squeeze(responseC(difpars.dif_dim(1),itime,jconds(icd,:)));
            rr2 = squeeze(responseC(difpars.dif_dim(2),itime,jconds(icd,:)));
            plot(rr1,rr2,'wo-');
        end
        if ~isfield(difpars,'axes') || isempty(difpars.axes)
            xlabel('dim1'); ylabel('dim2');
        else
            lx=xlabel(difpars.axes{1}); ly=ylabel(difpars.axes{2});
            set([lx ly],'interpreter','none');
        end
        hc=colorbar; axes(hc); ylabel('log10 counts');
        
        ha(4)=subplot(2,2,4); hi(4)=imagesc(ax,ay,log10(hpos(:,:)')); hold on;
        quiver(ax,ay,hvfd(:,:,1)'.*hmsk(:,:,1)',hvfd(:,:,2)'.*hmsk(:,:,2)','k');
        if ~isfield(difpars,'axes') || isempty(difpars.axes)
            xlabel('dim1'); ylabel('dim2');
        else
            lx=xlabel(difpars.axes{1}); ly=ylabel(difpars.axes{2});
            set([lx ly],'interpreter','none');
        end
        hc=colorbar; axes(hc); ylabel('log10 counts');
        
        matchc(hi(1:2),'abs');
        set(ha,'ydir','normal','plotbox',[1 1 1],...
            'xlim',[hx(1) hx(end)],'ylim',[hy(1) hy(end)]);
        scalesubplot(ha,0.07);
    end
else
    hf = [];
end

% Keep what you need
difdata.hist_response = hist_response;
difdata.hist_gradient = hist_gradient;
difdata.hist_gradientMask = hist_gradientMask;

return



