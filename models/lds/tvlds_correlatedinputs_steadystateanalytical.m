%{
---------------------------------------------------------------------------
This script computes analystical estimates of residual dynamics for a 
modularized "two-area" linear dynamical system, with an "input" area and a
"recorded" area. The "input" area provides inputs that are correlated across
time ("complex input"). We sample multiple models characterized by different
pairs of time-constants and rotation frequncies that are either assigned to
the input area or recorded area. We then obtain analytical estiamtes of the
residual dynamics assuming that only observations from the recorded area are
available. The analytical solutions help us to understand the 
inflationary effects of correlated input noise. See Supplementary Figure 4.
---------------------------------------------------------------------------
Author: Aniruddh Galgali (Nov 2020)

%}


clearvars -except DIRS
clc
close all
rseed = rng('default');
rng(rseed);


%% Estimation parameters

lags = [3]; % The iv lag provided to the analytical estimator (Supp Math Note B)

%% Simulation parameters

dt = 0.045; % Simulation step size (typically used for continuous time simulations)
obs_dim = 20; % number of observed dims for each area
latent_dim = 2; % number of latent dims for each area

eig_val_min = 0.4; % smallest eigenvalue that makes sense (this is tuned to a bin-size of 45 ms)
eig_val_max = 1.0 - 1e-4; % largest eigenvalue that is still stable. (We pick somethingnust below 1)

% Specify the range for the latent and observation noise (log-scale)
d_noise_latent_logmax = -5;
d_noise_latent_logmin = -7;



% Samping eigenvalues on a grid for both input and recorded areas.
% Specifies the sampling frequencey on the grid.
num_ev = 1000;
num_rot_f = 10;


% SCENARIOS
% {'1': input area has non-rotational dynamics and recorded area has 
%       rotational dynamics
%  '2': input area has rotational dynamics and recorded area has 
%       non-rotational dynamics
%  '3': both areas have decaying dynamics
% }

scenario_idx = 2;

if(scenario_idx == 1)
    
    input_dynamics = 'non-rotational';
    recurrent_dynamics = 'rotational';
    noise_eigen_basis_type = 'axis-aligned';
    dynamics_eigen_basis_type = '';
    
    tau_rot = 1;
    rot_f = [0.5 0.75 1.00 1.25 1.50];
    ev_samples = [0 linspace(eig_val_min, eig_val_max, num_ev)];
    n_models = length(ev_samples).*length(rot_f);

    vy_re_max = 0;
    vy_im_max = 1;
    vy_re_min = 0;
    vy_im_min = 1;
    
    ev_dyn = NaN(latent_dim,n_models);
    ev_noise = NaN(latent_dim,n_models); 
    V_dyn =  NaN(latent_dim,latent_dim,n_models);
    V_noise = repmat(eye(latent_dim),[1 1 n_models]);
    m_count = 1;
    
    for ir = 1:length(rot_f)
        for ie = 1: length(ev_samples)
            
            [A] = sample_single_rotational_dynamics_matrix(tau_rot, rot_f(ir),...
                vy_re_max, vy_re_min, vy_im_max, vy_im_min, dt, rseed);
            [vv,dd] = eig(A);
            ev_dyn(:,m_count) = diag(dd);
            V_dyn(:,:,m_count) = vv;
            ev_noise(:,m_count) = ev_samples(ie).*ones(latent_dim,1); 
            m_count = m_count + 1;
            
        end
        
    end
        
elseif(scenario_idx == 2)
    
    input_dynamics = 'rotational';
    recurrent_dynamics = 'non-rotational';
    dynamics_eigen_basis_type = 'axis-aligned';
    noise_eigen_basis_type = '';
    
    tau_rot = 1;
    rot_f = [0.5 0.75 1.00 1.25 1.50];
    ev_samples = [0 linspace(eig_val_min, eig_val_max, num_ev)];
    n_models = length(ev_samples).*length(rot_f);

    vy_re_max = 0;
    vy_im_max = 1;
    vy_re_min = 0;
    vy_im_min = 1;
    
    ev_dyn = NaN(latent_dim,n_models);
    ev_noise = NaN(latent_dim,n_models);
    V_noise =  NaN(latent_dim,latent_dim,n_models);
    V_dyn = repmat(eye(latent_dim),[1 1 n_models]);
    m_count = 1;
    for ir = 1:length(rot_f)
        for ie = 1: length(ev_samples)
            
            [A] = sample_single_rotational_dynamics_matrix(tau_rot, rot_f(ir),...
                vy_re_max, vy_re_min, vy_im_max, vy_im_min, dt, rseed);
            [vv,dd] = eig(A);
            ev_noise(:,m_count) = diag(dd);
            V_noise(:,:,m_count) = vv;
            ev_dyn(:,m_count) = ev_samples(ie).*ones(latent_dim,1);
            
            m_count = m_count + 1;
            
        end
        
    end
  
elseif(scenario_idx == 3)
    
    input_dynamics = 'non-rotational';
    recurrent_dynamics = 'non-rotational';
    dynamics_eigen_basis_type = 'axis-aligned';
    noise_eigen_basis_type = 'axis-aligned';
    tau_dyn_all = logspace(log10(0.05),log10(5),5);
    ev_samples = [0 linspace(eig_val_min, eig_val_max, num_ev)];

    n_models = length(tau_dyn_all).*length(ev_samples);
    ev_noise = NaN(latent_dim,n_models);
    ev_dyn = NaN(latent_dim,n_models);     
    V_dyn = repmat(eye(latent_dim),[1 1 n_models]);
    V_noise = repmat(eye(latent_dim),[1 1 n_models]);
    m_count = 1;
    
    for ir = 1: length(tau_dyn_all)
        for ie = 1: length(ev_samples)
            
            tau_dyn(1) = tau_dyn_all(ir);
            tau_dyn(2) = 0.4;
            ev_dyn(:,m_count) = convert_timeconstants_to_originaleigenvalues(tau_dyn, dt);
            ev_noise(:,m_count) = ev_samples(ie).*ones(latent_dim,1);

            m_count = m_count + 1;
        end
    end
 
end


% Sample latent and observation noise variances

d_noise_latent = NaN(latent_dim,n_models);


for ii = 1:n_models
  
    
    d_noise_latent(:,ii) = repmat((d_noise_latent_logmin + (d_noise_latent_logmax - ...
        d_noise_latent_logmin)*rand(1)),[latent_dim 1]);

end


%% Obtain the analytical solution for each model on the sampled grid for a given scenario

L = cell(length(lags),1);
for jj = 1:length(lags)
    out_labels = cell(n_models,1);
    for ii = 1:n_models

        xref = NaN(1,3*latent_dim+1);
        xref(1:latent_dim) = ev_dyn(:,ii)';
        xref(latent_dim+1 : 2*latent_dim) = ev_noise(:,ii)';
        xref(2*latent_dim+1 : 3*latent_dim) = 10.^ d_noise_latent(:,ii)';

        [out_labels{ii}] =  compute_analyticalsolutions_correlatedinputs(xref,latent_dim, dt,...
            V_dyn(:,:,ii), V_noise(:,:,ii),  lags(jj));

    end
    L{jj} = struct2table(cell2mat(out_labels));
end

%% PLOTTING THE ANALYTICAL ESTIMATRS ON THE GRID

col_reds = flipud(cbrewer('seq','Reds',15));

lag_idx_to_plot = 1;

if(scenario_idx == 1)
    
    figure;set(gcf,'Position',[1796  -95  1934  960]);
    
    eigval_eff_dim1 = L{lag_idx_to_plot}.eigval_eff_dim1;
    eigval_eff_dim2 = L{lag_idx_to_plot}.eigval_eff_dim2;
    
    eigval_noise_dim1 = L{lag_idx_to_plot}.eigval_noise_dim1;
    eigval_noise_dim2 = L{lag_idx_to_plot}.eigval_noise_dim2;
    rot_freq_dyn = L{lag_idx_to_plot}.rot_freq_dyn;
    
    
    ah = subplot(1,2,1); hold(ah); set(ah,'plotbox',[1 1 1]);
    unique_ev = unique(eigval_noise_dim1);
    ev_plot = [0 0.40 0.75 0.91 0.97 0.99];
    legend_labels = {};
    for ee = 1:length(ev_plot)
        [~,idx] = min(abs(unique_ev - ev_plot(ee)));
        idx_valid = find(abs(eigval_noise_dim1 - unique_ev(idx)) < 1e-4);
        plot(rot_freq_dyn(idx_valid),L{lag_idx_to_plot}.rot_freq_eff(idx_valid),'-','Color',col_reds(ee,:));
        legend_labels{ee} = sprintf('%s%.3f','ev=',ev_plot(ee));
        
    end
    legend(legend_labels);
    legend(legend_labels,'Location','NorthWest');
    xlabel('rot-freq-dyn')
    ylabel('rot-freq-eff')
    
    
    
    ah = subplot(1,2,2); hold(ah); set(ah,'plotbox',[1 1 1]);
    unique_rot_freq = unique(rot_freq_dyn );
    legend_labels = {};
    for ff = 1: length(unique_rot_freq)
        
        idx_valid = L{lag_idx_to_plot}.rot_freq_dyn == unique_rot_freq(ff);
        plot(ah,eigval_noise_dim1(idx_valid), eigval_eff_dim1(idx_valid),'-','Color',col_reds(ff,:));
        legend_labels{ff} = sprintf('%s%.2f','rotF=',unique_rot_freq(ff));
        
    end
    idx0 = eigval_noise_dim1 == ev_plot(1);
    plot(ah,[0.4 1],[uniquetol(eigval_eff_dim1(idx0),1e-3) uniquetol(eigval_eff_dim1(idx0),1e-3)], '--k','LineWidth',0.4);
    plot(ah,[uniquetol(eigval_eff_dim1(idx0),1e-3) uniquetol(eigval_eff_dim1(idx0),1e-3)],[0.93 1], '--k','LineWidth',0.4);
    hout = refline(1,0);
    set(hout,'Color',[0 0 0],'Parent',ah)
    set(ah,'xlim',[0.4 1]);
    set(ah,'ylim',[0.93 1]);
    legend(ah,legend_labels,'Location','NorthWest');
    xlabel(ah,'ev-noise')
    ylabel(ah,'ev-eff')
    
    ax_ticks = get(ah,'xtick');
    ay_ticks = get(ah,'ytick');
   
    ax1_pos = ah.Position; % position of first axes
    
    ax2 = axes('Position',ax1_pos,...
        'XAxisLocation','top',...
        'YAxisLocation','right',...
        'Color','none','plotbox',[1 1 1],'xlim',[0.4 1],'ylim',[0.93 1]);
    
    yticks_tau = abs(dt./log(ay_ticks));
    set(ax2,'YColor',[0 0 0])
    y_ticks_tau = arrayfun( @(x) sprintf('%s%.2f','-',x), yticks_tau,'uniformoutput', false);
    set(ax2,'yticklabel',y_ticks_tau);
    ylabel(ax2,'tau-eff')
    
    xticks_tau = -dt./log(ax_ticks);
    set(ax2,'XColor',[0 0 0])
    x_ticks_tau = arrayfun( @(x) sprintf('%s%.2f','-',x), xticks_tau,'uniformoutput', false);
    set(ax2,'xticklabel',x_ticks_tau);
    xlabel(ax2,'tau-noise')
    
elseif(scenario_idx == 2)

    figure;set(gcf,'Position',[1796  -95  1934  960]);
    
    eigval_eff_dim1 = L{lag_idx_to_plot}.eigval_eff_dim1;
    eigval_eff_dim2 = L{lag_idx_to_plot}.eigval_eff_dim2;
    
    eigval_dyn_dim1 = L{lag_idx_to_plot}.eigval_dyn_dim1;
    eigval_dyn_dim2 = L{lag_idx_to_plot}.eigval_dyn_dim2;
    rot_freq_inp = L{lag_idx_to_plot}.rot_freq_noise;
    
    ah = subplot(1,2,1); hold(ah); set(ah,'plotbox',[1 1 1]);
    unique_ev = unique(eigval_dyn_dim1);
    ev_plot = [0 0.40 0.75 0.91 0.97 0.99];
    legend_labels = {};
    for ee = 1:length(ev_plot)
       [~,idx] = min(abs(unique_ev - ev_plot(ee)));
       idx_valid = find(abs(eigval_dyn_dim1 - unique_ev(idx)) < 1e-4);
       plot(rot_freq_inp(idx_valid),L{lag_idx_to_plot}.rot_freq_eff(idx_valid),'-','Color',col_reds(ee,:));
       legend_labels{ee} = sprintf('%s%.2f','ev=',ev_plot(ee));
       
    end
    legend(legend_labels);
    legend(legend_labels,'Location','NorthWest');
    xlabel('rot-freq-noise')
    ylabel('rot-freq-eff')
    
    
    ah = subplot(1,2,2); hold(ah); set(ah,'plotbox',[1 1 1]);
    unique_rot_freq = unique(rot_freq_inp );
    
    legend_labels = {};
    for ff = 1: length(unique_rot_freq)
        
        idx_valid = L{lag_idx_to_plot}.rot_freq_noise == unique_rot_freq(ff);
        plot(eigval_dyn_dim1(idx_valid), eigval_eff_dim1(idx_valid),'-','Color',col_reds(ff,:));
        legend_labels{ff} = sprintf('%s%.2f','rotF=',unique_rot_freq(ff));

    end
    idx0 = eigval_dyn_dim1 == ev_plot(1);
    plot([0.4 1],[uniquetol(eigval_eff_dim1(idx0),1e-3) uniquetol(eigval_eff_dim1(idx0),1e-3)], '--k','LineWidth',0.4);
    plot([uniquetol(eigval_eff_dim1(idx0),1e-3) uniquetol(eigval_eff_dim1(idx0),1e-3)],[0.93 1], '--k','LineWidth',0.4);
    hout = refline(1,0);
    set(gca,'xlim',[0.4 1]);
    set(gca,'ylim',[0.93 1]);
    set(hout,'Color',[0 0 0])
    legend(legend_labels);
    legend(legend_labels,'Location','NorthWest');
    xlabel('ev-dyn')
    ylabel('ev-eff')
   
    ax_ticks = get(ah,'xtick');
    ay_ticks = get(ah,'ytick');
   
    ax1_pos = ah.Position; % position of first axes
    
    ax2 = axes('Position',ax1_pos,...
        'XAxisLocation','top',...
        'YAxisLocation','right',...
        'Color','none','plotbox',[1 1 1],'xlim',[0.4 1],'ylim',[0.93 1]);
    
    yticks_tau = abs(dt./log(ay_ticks));
    set(ax2,'YColor',[0 0 0])
    y_ticks_tau = arrayfun( @(x) sprintf('%s%.2f','-',x), yticks_tau,'uniformoutput', false);
    set(ax2,'yticklabel',y_ticks_tau);
    ylabel(ax2,'tau-eff')
    
    xticks_tau = -dt./log(ax_ticks);
    set(ax2,'XColor',[0 0 0])
    x_ticks_tau = arrayfun( @(x) sprintf('%s%.2f','-',x), xticks_tau,'uniformoutput', false);
    set(ax2,'xticklabel',x_ticks_tau);
    xlabel(ax2,'tau-noise')
    
elseif(scenario_idx == 3)
    
    figure;set(gcf,'Position',[2079         -85        1276         933]);

    eigval_eff_dim1 = L{lag_idx_to_plot}.eigval_eff_dim1;
    eigval_eff_dim2 = L{lag_idx_to_plot}.eigval_eff_dim2;
    
    eigval_noise_dim1 = L{lag_idx_to_plot}.eigval_noise_dim1;
    eigval_noise_dim2 = L{lag_idx_to_plot}.eigval_noise_dim2;
    
    eigval_dyn_dim1 = L{lag_idx_to_plot}.eigval_dyn_dim1;
    eigval_dyn_dim2 = L{lag_idx_to_plot}.eigval_dyn_dim2;
    
    unique_eigval_dyn = unique(eigval_dyn_dim1);
        
    ah = gca;hold(ah);set(ah,'plotbox',[1 1 1]);
    legend_labels = {};
    for dd = 1 : length(unique_eigval_dyn)
        
        idx_valid = find(eigval_dyn_dim1 == unique_eigval_dyn(dd));
        plot(eigval_noise_dim1(idx_valid(2:end)), eigval_eff_dim1(idx_valid(2:end)),'Color',col_reds(dd,:));
        legend_labels{dd} = sprintf('%s%.2f','evdyn=',uniquetol(eigval_dyn_dim1(idx_valid),1e-3));

    end
    set(gca,'xlim',[0.4 1]);
    set(gca,'ylim',[0.4 1]);
    legend(legend_labels,'Location','SouthEast');
    xlabel('ev-noise')
    ylabel('ev-eff')
    for dd = 1 : length(unique_eigval_dyn)
        
        idx_valid = find(eigval_dyn_dim1 == unique_eigval_dyn(dd));
        plot(eigval_noise_dim1(idx_valid(2)), eigval_eff_dim1(idx_valid(1)),'o',...
            'MarkerEdgeColor',col_reds(dd,:),'MarkerFaceColor',col_reds(dd,:));
        
    end
    
    ah = gca;
 
    ax_ticks = get(ah,'xtick');
    ay_ticks = get(ah,'ytick');
   
    ax1_pos = ah.Position; % position of first axes
    
    ax2 = axes('Position',ax1_pos,...
        'XAxisLocation','top',...
        'YAxisLocation','right',...
        'Color','none','plotbox',[1 1 1],'xlim',[0.4 1],'ylim',[0.4 1]);
    
    yticks_tau = abs(dt./log(ay_ticks));
    set(ax2,'YColor',[0 0 0])
    y_ticks_tau = arrayfun( @(x) sprintf('%s%.2f','-',x), yticks_tau,'uniformoutput', false);
    set(ax2,'yticklabel',y_ticks_tau);
    ylabel(ax2,'tau-eff')
    
    xticks_tau = -dt./log(ax_ticks);
    set(ax2,'XColor',[0 0 0])
    x_ticks_tau = arrayfun( @(x) sprintf('%s%.2f','-',x), xticks_tau,'uniformoutput', false);
    set(ax2,'xticklabel',x_ticks_tau);
    xlabel(ax2,'tau-noise')

end
