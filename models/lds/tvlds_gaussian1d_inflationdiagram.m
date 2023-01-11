%{
---------------------------------------------------------------------------
This script simulates a time-variying dynamical system (gaussian ibs.) with
a single dimensional latent state to understand the inflationary effects of
correlated input noise and contrast it to switches in the time course of 
latent noise. See Supplementary Figure 5.
---------------------------------------------------------------------------
Author: Aniruddh Galgali (Dec 2020)
%}

clearvars -except DIRS
clc
rseed = rng('default');
rng(rseed);
save_fig_format = '-pdf';
save_fig_path = './figures/simulations'; % create this path if it doesnt exist
%%
options.rseed = rseed;

% For this analysis we only consider a discrete dynamical process
options.dynamics_model = 'discrete';

options.T = 700; % Number of time points
options.lat_dim = 2;
options.num_trials = 1000;
options.dt = 0.045; % simulation step size


%% Generate dynamical parameters


% {'1': steady-state, no input correlations, fixed noise
%  '2': no input correlations, fixed noise -> input correlations, fixed noise
%  '3': no input correlations, switching noise
% }

scenarioID = 3;
time_of_switch = 175;

if(scenarioID == 1)

    % Setting model parameters  
    tau_recorded = 0.3;
    tau_input = 0;
    
    ev_recorded = convert_timeconstants_to_originaleigenvalues(tau_recorded,options.dt);
    latent_noise= 1e-6; % 1e-8 and 1e-11 for correlated input
    ev_input = convert_timeconstants_to_originaleigenvalues(tau_input,options.dt);

    % Generate parameter matrices and compute steady-state covariance
    A_augmented = [ev_recorded 1; 0 ev_input];
    Q_augmented = [1e-30 0 ; 0 latent_noise];
    
    Q0 = dlyap(A_augmented, Q_augmented);
    A_t = repmat(A_augmented,[1 1 options.T]);
    Q_t = repmat(Q_augmented,[1 1 options.T]);
    
    time_idxs_to_plot = [options.T - 50];
    
elseif(scenarioID == 2)
    
    % Setting model parameters  
    tau_recorded = 0.3;
    tau_input_pre = 0;
    tau_input_post = 0.2;
    
    ev_recorded = convert_timeconstants_to_originaleigenvalues([tau_recorded],options.dt);    
    ev_input_pre = convert_timeconstants_to_originaleigenvalues(tau_input_pre,options.dt);
    ev_input_post = convert_timeconstants_to_originaleigenvalues(tau_input_post,options.dt);

    % Generate parameter matrices and compute steady-state covariance

    A_augmented_pre = [ev_recorded 1;0 ev_input_pre];
    A_augmented_post = [ev_recorded 1;0 ev_input_post];
    
    latent_noise= 1e-6;
    Q_augmented = [1e-30 0 ; 0 latent_noise];
  
    Q0 = dlyap(A_augmented_pre, Q_augmented);
    Q_t = repmat(Q_augmented,[1 1 options.T]);
    
    A_t = cat(3, repmat(A_augmented_pre,[1 1 time_of_switch]),...
        repmat(A_augmented_post,[1 1 options.T - time_of_switch]));
    
    time_idxs_to_plot = [time_of_switch + 5 options.T - 50];
    
elseif(scenarioID == 3)
    
    % Setting model parameters  
    tau_recorded = 0.3;
    tau_input = 0;
    
    ev_recorded = convert_timeconstants_to_originaleigenvalues(tau_recorded, options.dt);
    ev_input = convert_timeconstants_to_originaleigenvalues(tau_input, options.dt);

    latent_noise_pre = 1e-6; % 1e-8 and 1e-11 for correlated input
    latent_noise_post = 1e-5;
    Q_augmented_pre = [1e-30 0 ; 0 latent_noise_pre];
    Q_augmented_post = [1e-30 0 ; 0 latent_noise_post];

    % Generate parameter matrices and compute steady-state covariance
    A_augmented = [ev_recorded 1; 0 ev_input];
    
    Q0 = dlyap(A_augmented, Q_augmented_pre);
    A_t = repmat(A_augmented,[1 1 options.T]);
    
    Q_t = cat(3, repmat(Q_augmented_pre,[1 1 time_of_switch]),...
    repmat(Q_augmented_post,[1 1 options.T - time_of_switch]));
    
    time_idxs_to_plot = [time_of_switch + 5 options.T - 50];
    
else
    
    error('invalid scenario')
    
end
options.A = A_t;
options.Q = Q_t;
options.sigma0 = Q0;

%% Generate latent states

% Since we are only ever using the latent states directly for this
% analysis, I don't simulate noisy observations. Therfore, I directly use 
% the function that generates latent states

inputs = zeros(options.lat_dim,options.T,options.num_trials);
options.x0 = zeros(options.lat_dim, options.num_trials);
for tt= 1:options.T

    inputs(:,tt,:) = chol(options.Q(:,:,tt))'*randn(options.lat_dim,options.num_trials);

end

[latent_states] = generateLatentStates(options.A, options.x0, options.sigma0,...
     inputs, options.dynamics_model, options.dt, options.rseed);
propagated_state = latent_states(1,2:end,:) - latent_states(2,1:end-1,:);


%% Plot scatter plots of latent vs propagated states to demonstrate inflation

plotpars.marker = 'o';
plotpars.linewidth = 0.1;
plotpars.marker_edge_col = [0 0 0];
plotpars.marker_face_col = [0.5 0.5 0.5];
plotpars.marker_edge_alpha = 0.5;
plotpars.marker_face_alpha = 1;
plotpars.marker_size = 12;

    
figure;
fit_line = fittype('a*x+b');
for tt = 1: length(time_idxs_to_plot)

    ah = subplot(length(time_idxs_to_plot),3, (tt-1)* 3 + 1);hold(ah);set(ah,'plotbox',[1 1 1]);
    set(ah,'xlim',[-0.025 0.025],'xtick',[-0.02 -0.01 0 0.01 0.02],'ylim',[-0.025 0.025],'ytick',[-0.02 -0.01 0 0.01 0.02]);    
    scatter(ah, squeeze(latent_states(1,time_idxs_to_plot(tt),:)),squeeze(propagated_state(1,time_idxs_to_plot(tt),:)),...
        plotpars.marker_size, 'LineWidth', plotpars.linewidth,'MarkerEdgeColor',plotpars.marker_edge_col,...
        'MarkerFaceColor', plotpars.marker_face_col,'MarkerEdgeAlpha',plotpars.marker_edge_alpha,...
        'MarkerFaceAlpha',plotpars.marker_face_alpha); 
    xlabel(sprintf('%s%d','x(t) : t = ', time_idxs_to_plot(tt)),'fontsize',8);
    ylabel(sprintf('%s%d','A(t)*x(t) : t = ', time_idxs_to_plot(tt)),'fontsize',8);
    xlims = get(gca,'xlim');
    line([xlims(1) xlims(end)],[xlims(1) xlims(end)],'Color','red')
    fm = fit(squeeze(latent_states(1,time_idxs_to_plot(tt),:)),squeeze(propagated_state(1,time_idxs_to_plot(tt),:)),fit_line);
    xFit = linspace(xlims(1), xlims(end), 500);
    yFit = fm.a*xFit + fm.b;
    plot(xFit, yFit, '-', 'Color',[0.5 0.5 0.5]); % Plot fitted line.
    title(sprintf('%s%.6f','reg-coef = ',fm.a)); 

    ah = subplot(length(time_idxs_to_plot),3,(tt-1)* 3 + 2);hold(ah);set(ah,'plotbox',[1 1 1]);
    set(ah,'xlim',[-0.025 0.025],'xtick',[-0.02 -0.01 0 0.01 0.02],'ylim',[-0.025 0.025],'ytick',[-0.02 -0.01 0 0.01 0.02]); 
    scatter(ah, squeeze(latent_states(1,time_idxs_to_plot(tt),:)),squeeze(latent_states(2,time_idxs_to_plot(tt),:)),...
    plotpars.marker_size, 'LineWidth', plotpars.linewidth,'MarkerEdgeColor',plotpars.marker_edge_col,...
    'MarkerFaceColor', plotpars.marker_face_col,'MarkerEdgeAlpha',plotpars.marker_edge_alpha,...
    'MarkerFaceAlpha',plotpars.marker_face_alpha); 
    xlabel(sprintf('%s%d','x(t) : t = ', time_idxs_to_plot(tt)),'fontsize',8);
    ylabel(sprintf('%s%d','eps(t) : t = ', time_idxs_to_plot(tt)),'fontsize',8);
    ylims = get(gca,'ylim');
    line([ylims(1) ylims(end)],[ylims(1) ylims(end)],'Color','red');
    fm = fit(squeeze(latent_states(1,time_idxs_to_plot(tt),:)),squeeze(latent_states(2,time_idxs_to_plot(tt),:)),fit_line);
    xFit = linspace(xlims(1), xlims(end), 500);
    yFit = fm.a*xFit + fm.b;
    plot(xFit, yFit, 'Color',[0.5 0.5 0.5]); % Plot fitted line.
    title(sprintf('%s%.6f','reg-coef = ',fm.a)); 

    ah = subplot(length(time_idxs_to_plot),3,(tt-1)* 3 + 3);hold(ah);set(ah,'plotbox',[1 1 1]);
    set(ah,'xlim',[-0.025 0.025],'xtick',[-0.02 -0.01 0 0.01 0.02],'ylim',[-0.025 0.025],'ytick',[-0.02 -0.01 0 0.01 0.02]); 
    scatter(ah, squeeze(latent_states(1,time_idxs_to_plot(tt),:)),squeeze(latent_states(1,time_idxs_to_plot(tt) + 1,:)),...
     plotpars.marker_size, 'LineWidth', plotpars.linewidth,'MarkerEdgeColor',plotpars.marker_edge_col,...
    'MarkerFaceColor', plotpars.marker_face_col,'MarkerEdgeAlpha',plotpars.marker_edge_alpha,...
    'MarkerFaceAlpha',plotpars.marker_face_alpha);     
    xlabel(sprintf('%s%d','x(t) : t = ', time_idxs_to_plot(tt)),'fontsize',8);
    ylabel(sprintf('%s%d','x(t+1) : t = ', time_idxs_to_plot(tt)),'fontsize',8);
    xlims = get(gca,'xlim');
    line([xlims(1) xlims(end)],[xlims(1) xlims(end)],'Color','red')
    fm = fit(squeeze(latent_states(1,time_idxs_to_plot(tt),:)),squeeze(latent_states(1,time_idxs_to_plot(tt)+1,:)),fit_line);
    xFit = linspace(xlims(1), xlims(end), 500);
    yFit = fm.a*xFit + fm.b;
    plot(xFit, yFit, 'Color',[0.5 0.5 0.5]); % Plot fitted line.
    title(sprintf('%s%.6f','reg-coef = ',fm.a));
end

%% Save figures

set(gcf,'Position',[322   420   932   619]);
fig_name = 'EVInflationdiagram_1Ddynamics_supplementfigure';
export_fig(fullfile(save_fig_path,fig_name),'-painters','-transparent',save_fig_format)
