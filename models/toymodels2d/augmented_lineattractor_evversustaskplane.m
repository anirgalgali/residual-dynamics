%{ This script simulates an augmented line attractor model, which
%  contains 2 additional dimensions. Both of these additional dimensions
%  are associated with quickly decaying reurrent dynamics but driven by
%  strong sinusoidal inputs that re condition-independent (see Figure 5).
%
%}

clearvars -except DIRS
clc
close all

%%
data_load_path = fullfile(DIRS.analysis,'/simulations/toymodels/');
file_name = 'integration_models_newAnalysisPipeline_09-Dec-2018.mat';
load(fullfile(data_load_path,file_name));
save_fig_path = '';% specify a valid path for saving figures.
save_file_path = fullfile(DIRS.analysis,'/analyses/');
do_save_fig = false;

%% Setting high-level simulation parameters

rseed =  model_analyses.simPars.model2.model_obj.RAND_SEED; % Taken from Toy DynamicsModel.m 
rng(rseed);

% Latent dims
n_choice_dep_dims = 2;
n_choice_idp_dims = 2;
n_total_latent_dims = n_choice_dep_dims + n_choice_idp_dims;

% Simulation parameters
n_trials = model_analyses.simPars.model2.model_obj.num_trajectories;
n_times = model_analyses.simPars.model2.model_obj.sim_time;
dt = model_analyses.simPars.model2.model_obj.TIME_STEP;

% parameters determining inital condition of latent states

x_latent_noise_var = unique(diag(model_analyses.simPars.model2.model_obj.noise_var));

% parameters determining choice-dependent input time-course
choice_dep_input_magnitude = []; % empty if the inputs are loaded from disk (default)
choice_dep_input_time_onset = 0.0;
choice_dep_input_time_offset = 1.0;
choice_dep_rot_freq = [];

% parameters determining choice-independent rotational input time-course
choice_idp_input_magnitude = 20;
choice_idp_input_time_onset = 0;
choice_idp_input_time_offset = 1.0;
choice_idp_rot_freq = 0.8; % units of Hz

% parameters determmining dynamics of the latent states

pars_nonNormal_la.A = model_analyses.simPars.model2.pars.A;
lambda_decay_1 = -7.5;
lambda_decay_2 = -10;


% parameter determining overlap between condition-independent and
% condition-dependent modes
angle_between_planes = 90;

% parameters determinng observation model and binning of time-series
obs_dim = n_total_latent_dims;
do_orth_obs_mat = true;
obs_noise_var_knob = 1e-6;
obs_noise_var = diag(obs_noise_var_knob.*rand(obs_dim,1));
bin_size = [45];

% choice labels and inputs of stored simulation variables
choice_labels = model_analyses.data.states.model2.choice_labels;
model_optimized_inputs = model_analyses.data.inputs.model2.mean;

%% Generating inputs to model

times = 0:dt:(n_times-1)*dt;

% GENERATING CONDITION-DEPENDENT INPUTS

% Genearting condition dependent input time masks
choice_dep_input_mask = false(1,n_times-1);
choice_dep_input_times_on = times(1:end-1) >= choice_dep_input_time_onset & ...
    times(1:end-1) <= choice_dep_input_time_offset;
choice_dep_input_mask(choice_dep_input_times_on) = true;
% 

% choice_dep_input_magnitude_t = repmat(choice_dep_input_magnitude.*ones(n_choice_dep_dims,1),[1 n_times - 1]);
choice_dep_input_t = NaN(n_choice_dep_dims, n_times-1, n_trials);
for kk = 1: n_trials
    
    if(choice_labels(kk) == -1)
        cidx = 1;
    elseif(choice_labels(kk) == +1)
        cidx = 2;
    end
        
    for tt = 1 : n_times - 1
        choice_dep_input_t(:,tt,kk) = model_optimized_inputs(:,tt,cidx);
    end
end
% choice_dep_input_t = inputs_model2;
choice_dep_input_t = choice_dep_input_mask.* choice_dep_input_t;

% GENERATING CONDITION-INDEPENDENT INPUTS

% Genearting condition independent input time masks
choice_idp_input_mask = false(1,n_times);
choice_idp_input_times_on = times >= choice_idp_input_time_onset & ...
    times <= choice_idp_input_time_offset;
choice_idp_input_mask(choice_idp_input_times_on) = true;

% Genearting time course of condition-independent inputs

choice_idp_input_t = NaN(n_choice_idp_dims, n_times, n_trials);
for tt = 1 : n_times
    
    choice_idp_input_t(:,tt,:) = repmat(choice_idp_input_magnitude.*([cos(2*pi*choice_idp_rot_freq * times(tt));...
                                                                      sin(2*pi*choice_idp_rot_freq * times(tt))]),[1 n_trials]);      
end
choice_idp_input_t = choice_idp_input_mask.*choice_idp_input_t;

% Final inputs to dynamics
inputs_to_dynamics = cat(1,choice_dep_input_t,choice_idp_input_t(:,1:end-1,:));


%% Generating latent dynamics


latent_states = NaN(n_total_latent_dims, n_times, n_trials);
latent_states_noiseless = NaN(size(latent_states));
xdots_sim = NaN(size(latent_states));
latent_noise_cov = diag(x_latent_noise_var .* ones(n_total_latent_dims,1));

% Constructing the full dynamics matrix

[v_eig_la, d_eig_la] = eig(pars_nonNormal_la.A);
V_full = blkdiag(v_eig_la, eye(n_choice_idp_dims));
A_full = V_full * diag([diag(d_eig_la); lambda_decay_1; lambda_decay_2]) * inv(V_full);

% Determining steady-state initial variance for the augmented dimensions

rel_dims = n_choice_dep_dims + 1 : n_choice_dep_dims + n_choice_idp_dims;
[C_ss0] = compute_augmentedLDS_correlatedInputNoise(A_full(rel_dims,rel_dims), dt,ones(length(rel_dims),1), (1./dt.^2).*latent_noise_cov(rel_dims, rel_dims));

% Definining initial conditions

x0_original_dims = squeeze(model_analyses.data.states.model2.response(:, 1 , :));
x0_augmented_dims = chol(C_ss0(1 : floor(size(C_ss0,1)/2),1:floor(size(C_ss0,1)/2)))*randn(n_choice_idp_dims, n_trials);
x0 = cat(1,x0_original_dims,x0_augmented_dims);

% Generating latent state trajectories
for tt = 1: n_times
    
    if(tt == 1)
         latent_states(:,tt,:) = x0;
         latent_states_noiseless(:,tt,:) = x0;
    end
    
    if(tt < n_times) 
        x_dot_instantaneous = A_full * squeeze(latent_states(:,tt,:));
        xdots_sim(:,tt,:) = x_dot_instantaneous;
        latent_states_noiseless(:,tt+1,:) = squeeze(latent_states(:,tt,:)) + dt.*(x_dot_instantaneous + squeeze(inputs_to_dynamics(:,tt,:)));
        latent_states(:,tt+1,:)  = squeeze(latent_states_noiseless(:,tt+1,:)) + chol(latent_noise_cov)'*randn(n_total_latent_dims, n_trials);
    end
end
%% Plotting latent state trajectories

n_trajs_to_plot = 500;
idx_trials_to_plot = randperm(n_trials,n_trajs_to_plot);

choice_labels_of_plotted_trials = choice_labels(idx_trials_to_plot);
unique_choice_labels = unique(choice_labels_of_plotted_trials);
trajcols = {'-b','-r'};

figure; ah1 = subplot(1,2,1);hold(ah1);set(ah1,'plotbox',[1 1 1])
for ichoice = 1: length(unique_choice_labels)
    latent_states_to_plot = latent_states(:,:,idx_trials_to_plot);
    plot(ah1, squeeze(latent_states_to_plot(1,:,choice_labels_of_plotted_trials == unique_choice_labels(ichoice))),...
          squeeze(latent_states_to_plot(2,:,choice_labels_of_plotted_trials == unique_choice_labels(ichoice))),trajcols{ichoice})
end
xlabel(ah1,'x_{1}');
ylabel(ah1,'x_{2}');
axis tight
ah2 = subplot(1,2,2);hold(ah2);set(ah2,'plotbox',[1 1 1])
for ichoice = 1: length(unique_choice_labels)
    latent_states_to_plot = latent_states(:,:,idx_trials_to_plot);
    plot(ah2, squeeze(latent_states_to_plot(3,:,choice_labels_of_plotted_trials == unique_choice_labels(ichoice))),...
          squeeze(latent_states_to_plot(4,:,choice_labels_of_plotted_trials == unique_choice_labels(ichoice))),trajcols{ichoice})
end
xlabel(ah2,'x_{3}');
ylabel(ah2,'x_{4}');
axis tight

if(do_save_fig)
    export_fig(fullfile(save_fig_path, sprintf('%s%d%s%.1f%s%.1f','latentstatetrajs_augmentedlineatt_planeangle=',...
        angle_between_planes,'_lambda_decay=',lambda_decay_1,'_',lambda_decay_2)),'-painters','-transparent','-pdf')
end

figure;
coltrajs_c1 = flipud(cbrewer('seq','GnBu',6));
coltrajs_c2 = flipud(cbrewer('seq','OrRd',6));
ah3 = subplot(1,2,1);hold(ah3);set(ah3,'plotbox',[1 1 1]);
for idim = 1:n_total_latent_dims
    plot(ah3, times, squeeze(var(latent_states(idim,:,choice_labels == -1),0,3)),'Color',coltrajs_c1(idim,:))
end

ah4 = subplot(1,2,2);hold(ah4);set(ah4,'plotbox',[1 1 1]);
for idim = 1:n_total_latent_dims
    plot(ah4, times, squeeze(var(latent_states(idim,:,choice_labels == +1),0,3)),'Color',coltrajs_c2(idim,:))
end

%% Generating observations

clear model_analyses % Clearing up unnecessary variables

% Computing observation matrices 

cos_theta = @(theta) diag(cos(deg2rad(90 - theta)) .* ones(n_total_latent_dims/2, 1));
sin_theta = @(theta) diag(sin(deg2rad(90 - theta)) .* ones(n_total_latent_dims/2, 1));
rot_matrix = [cos_theta(angle_between_planes) -1.*sin_theta(angle_between_planes) ; ...
              sin_theta(angle_between_planes) cos_theta(angle_between_planes)]; 
          
if(obs_dim == n_total_latent_dims)
    C = eye(obs_dim);
else
    if(do_orth_obs_mat)
        C = RandOrthMat(obs_dim);
    else
        C = randn(obs_dim);
    end
    C = C(:,1:n_total_latent_dims);
end

Crot = C*rot_matrix;
C = [C(:,1:2) Crot(:,3:4)];

% Generating observations
observed_states = NaN(obs_dim, n_times, n_trials);
for tt = 1: n_times
    observed_states(:,tt,:) = C*squeeze(latent_states(:,tt,:)) + chol(obs_noise_var)*randn(obs_dim, n_trials);
end

%% Storing all relevant simulation variables

sim_opts.rseed = rseed;

sim_opts.dim_pars.n_choice_dep_dims = n_choice_dep_dims;
sim_opts.dim_pars.n_choice_idp_dims = n_choice_idp_dims;
sim_opts.dim_pars.n_total_latent_dims = n_total_latent_dims;
sim_opts.dim_pars.obs_dim = obs_dim;

sim_opts.expt_pars.n_trials = n_trials;
sim_opts.expt_pars.n_times = n_times;
sim_opts.expt_pars.sim_dt = dt;
sim_opts.expt_pars.bin_size = bin_size;

sim_opts.input_pars.choice_dep.magnitude = choice_dep_input_magnitude;
sim_opts.input_pars.choice_dep.time_onset = choice_dep_input_time_onset;
sim_opts.input_pars.choice_dep.time_offset = choice_dep_input_time_offset;
sim_opts.input_pars.choice_dep.input_mask = choice_dep_input_mask;
sim_opts.input_pars.choice_dep.rot_freq = choice_dep_rot_freq;

sim_opts.input_pars.choice_idp.magnitude = choice_idp_input_magnitude;
sim_opts.input_pars.choice_idp.time_onset = choice_idp_input_time_onset;
sim_opts.input_pars.choice_idp.time_offset = choice_idp_input_time_offset;
sim_opts.input_pars.choice_idp.input_mask = choice_idp_input_mask;
sim_opts.input_pars.choice_idp.rot_freq = choice_idp_rot_freq;

sim_opts.latent_model.model_pars.A = A_full;
sim_opts.latent_model.model_pars.Q = latent_noise_cov;

sim_opts.obs_model.model_pars.C = C;
sim_opts.obs_model.model_pars.R = obs_noise_var;
sim_opts.obs_model.do_orth_obs_mat = do_orth_obs_mat;
sim_opts.obs_model.angle_between_planes = angle_between_planes;

sim_opts.trial_info.x_ics = x0;
sim_opts.trial_info.choice_labels = choice_labels;

 %% Binning data and computing residuals

clear data_
data_.response = observed_states;
trial_labels = zeros(n_trials,1);
trial_labels(choice_labels == -1) = 1;
trial_labels(choice_labels == +1) = 2;
data_.task_variable.targ_dir = trial_labels;
[~,~,data_.task_index.targ_dir] = unique(trial_labels);
data_.time = times;
data_.time_rel = times;
data_.time_iev = ones(1,length(times));
data_ = binSpikeCounts(data_, bin_size,1./dt, false, false);
residual_condition_vars{1} = {'targ_dir'};
[data_res] = compute_residual_trajectories({data_},residual_condition_vars);

[sim_opts.latent_model.binned_model_pars.A] = ...
    compute_binned_groundtruthdynamics(sim_opts.latent_model.model_pars.A,...
    sim_opts.expt_pars.sim_dt, sim_opts.expt_pars.bin_size);

%% Computing task-relevant planes and their projections

alignments = {'stim'};

jPC_pars.n_jpc_planes = 2;
jPC_pars.do_pre_smooth = false;
jPC_pars.time_win_analyze = {[0 1]}; 
jPC_pars.alignments = alignments;

choice_time_pars.pref_idx = 1;
choice_time_pars.anti_idx = 2;
choice_time_pars.condition_vars = {'targ_dir'};
choice_time_pars.alignments = alignments;
[analysis_condavg] = extract_conditionaverage_planes_v2(data_, trial_labels, jPC_pars, choice_time_pars, alignments);
%%
scale_ = 250;
gridx = -100:5:100;
gridy = -100:5:100;
flow_pars.start_X = 1;
flow_pars.start_Y = 1;
flow_pars.end_X = length(gridx) ; 
flow_pars.end_Y = length(gridy);
flow_pars.step_X = 5;
flow_pars.step_Y = 5;
flow_pars.xLim =  [-110 110];
flow_pars.yLim =  [-110 110];
flow_pars.xTick =  [-100:25:100];
flow_pars.yTick=  [-100:25:100];
flow_pars.xTickLabel =  [-100:25:100];
flow_pars.yTickLabel =  [-100:25:100];
flow_pars.figPosition = [130 69 635 636];
flow_pars.Color = [0.75 0.75 0.75];
flow_pars.dotColor = [0 0 0];
flow_pars.dotSize = 4;
flow_pars.lineWidth  = 1;
flow_pars.LineLength = 0.02;
flow_pars.arrowShape = [1 0.6 scale_*0.004 scale_*0.06 scale_*0.01 scale_*0.025 scale_*0.000 scale_*0.005];
flow_pars.plot_locs = false;
for ialign = 1: length(alignments)
    h = figure; 
    for ii = 1: length(analysis_condavg.(alignments{ialign}).projections)
       
        ah = subplot(2,2,ii); hold(ah); set(ah,'plotbox',[1 1 1])
        

        proj_plane = real(analysis_condavg.(alignments{ialign}).directions{ii});

        [locs, flow] = compute_flowfield_arbitraryplane(sim_opts.latent_model.binned_model_pars.A,sim_opts.obs_model.model_pars.C,...
            proj_plane, gridx, gridy, 'discrete_time',sim_opts.expt_pars.bin_size*sim_opts.expt_pars.sim_dt);

        [h,ah] = plotFlowField_v2(locs,flow,flow_pars,h,ah); 
        
        plot(ah,squeeze(analysis_condavg.(alignments{ialign}).projections{ii}(1,:,1)),...
            squeeze(analysis_condavg.(alignments{ialign}).projections{ii}(2,:,1)),'-b');
       
        plot(ah,squeeze(analysis_condavg.(alignments{ialign}).projections{ii}(1,:,1)),...
            squeeze(analysis_condavg.(alignments{ialign}).projections{ii}(2,:,1)),'o',...
            'markerfacecolor',[1 1 1],'markeredgecolor',[0 0 1]);
        
        plot(ah,squeeze(analysis_condavg.(alignments{ialign}).projections{ii}(1,:,2)),...
            squeeze(analysis_condavg.(alignments{ialign}).projections{ii}(2,:,2)),'-r');
        
        plot(ah,squeeze(analysis_condavg.(alignments{ialign}).projections{ii}(1,:,2)),...
            squeeze(analysis_condavg.(alignments{ialign}).projections{ii}(2,:,2)),'o',...
            'markerfacecolor',[1 1 1],'markeredgecolor',[1 0 0]);
        
        xlims = get(ah,'xlim');
        ylims = get(ah,'ylim');
        
        if(xlims(2) - xlims(1) > ylims(2) - ylims(1))
            set(ah,'xlim',xlims,'ylim',xlims);
        elseif(xlims(2) - xlims(1) < ylims(2) - ylims(1))
            set(ah,'xlim',ylims,'ylim',ylims);
        end
%            
        
        title(ah,strcat(analysis_condavg.(alignments{ialign}).direction_labels{ii},': R_{2} = ', num2str(analysis_condavg.(alignments{ialign}).variance_explained(ii)),...
            ', rotf = ',num2str(unique(abs(analysis_condavg.(alignments{ialign}).rot_evals{ii}/(2*pi))))));
%         
    end
end
if(do_save_fig)
    export_fig(fullfile(save_fig_path, sprintf('%s%d%s%.1f%s%.1f','condavg2Dtaskplanes_augmentedlineatt_planeangle=',...
        angle_between_planes,'_lambda_decay=',lambda_decay_1,'_',lambda_decay_2)),'-painters','-transparent','-pdf')
end
%% Running the residual dynamics pipeline

opts_analysis.hankel_order = 5;
opts_analysis.hankel_cv.hankel_ranks = [1 2 3 4 5];
opts_analysis.hankel_cv.n_cv = 20;
opts_analysis.hankel_cv.criterion = 'min';
opts_analysis.hankel_cv.type= 'combined';
opts_analysis.lag_cv.n_cv = 1;
opts_analysis.lag_cv.n_folds = 10;
opts_analysis.lag_cv.grid_type = 'uniform';
opts_analysis.lag_cv.do_display = true;

switch opts_analysis.lag_cv.grid_type
    
    case 'uniform'
        opts_analysis.lag_cv.iv_lag = [1 2 3 4 5];
        opts_analysis.lag_cv.sub_dim = [1 2 3 4];
        
    case 'random'
        opts_analysis.lag_cv.numPars = 20;
        opts_analysis.lag_cv.max_lag = 10;
        opts_analysis.lag_cv.max_dim = 6;
             
end
opts_analysis.doFinalFit = true;
opts_analysis.doSmoothCV = true;
opts_analysis.final_fit.lag_cv_metric = 'mse';
opts_analysis.final_fit.use_common_lagdim_across_conditions = true;
if(~opts_analysis.doSmoothCV)
    opts_analysis.final_fit.smoothalpha = {1e4};
end
opts_analysis.smooth_cv.alpha_all = [1e-3 1e-1 1e2 1e4 1e5 1e6 1e8];
opts_analysis.smooth_cv.n_cv = 5;
opts_analysis.smooth_cv.n_folds = 5;
opts_analysis.smooth_cv.combine_align = true;
opts_analysis.final_fit.do_cross_validate_noise_model = false;

[analysis_resdyn] = run_pipeline(data_res{1}.response, data_res{1}.time, data_res{1}.time_rel, data_res{1}.time_iev, trial_labels, opts_analysis, rseed);


%% Plotting cross-validation results

clear plotpars
col_trains = flipud(cbrewer('seq','OrRd',10));
col_tests = flipud(cbrewer('seq','Greys',10));
plotpars.err_type = 'pvar';
plotpars.x_type = 'dim';
plotpars.err_scale = 1.0;
plotpars.err_metric = 'sem';
plotpars.linewidth = 1.0;
plotpars.bar_length = 0.5;
plotpars.bar_width = 1.0;
plotpars.do_plot_train = false;
plotpars.markersize = 4;

switch plotpars.x_type
    case 'dim'
        plotpars.col_test = col_tests(1:length(analysis_resdyn.fit_pars.lag_cv.iv_lag),:);
        plotpars.col_train = col_trains(1:length(analysis_resdyn.fit_pars.lag_cv.iv_lag),:);
    case 'lag'
        plotpars.col_test = col_tests(1:length(analysis_resdyn.fit_pars.lag_cv.sub_dim),:);
        plotpars.col_train = col_trains(1:length(analysis_resdyn.fit_pars.lag_cv.sub_dim),:);
end

h = plotCV_lagDim(analysis_resdyn.cv_lagdim, analysis_resdyn.cv_lagdim_hyp, plotpars);
ax = get(h,'Children');
ax = flipud(ax);
for iax = 1: length(ax)
    title(ax(iax),sprintf('%s%d%s%d','optLag = ',analysis_resdyn.cv_lagdim_opt{iax}.mse_opt.lags(1),...
        ', optDim =',analysis_resdyn.cv_lagdim_opt{iax}.mse_opt.dims(1)),'FontSize',8);
end
set(gcf,'Position',[566 -318 973 1131]);


%% Plotting the final fit results

figure; set(gcf, 'Position', [1469 44 1594 609]);
ah = subplot(1, 3, 1);set(ah,'plotbox',[1 1 1]);hold(ah)
 

clear plotpars
clear data_to_plot

time_abs = analysis_resdyn.result.final_model{1}.time;
col_lines = flipud(cbrewer('seq','GnBu',20));
plotpars.col_lines = col_lines(round(linspace(1, size(col_lines,1) - 10, size(analysis_resdyn.result.final_model{1}.A,1))),:);
plotpars.do_plot_markers = true;
plotpars.markersize = 4;
plotpars.linewidth = 2;
plotpars.ylim = [0 1.5];
plotpars.ytick = [0 : 0.2 : 1.5];
plotpars.yticklabel = [0 : 0.2 : 1.5];
plotpars.yLabel = '| eigen-value | (a.u)';
plotpars.xLabel = 'time (s)';
plotpars.do_plot_stability = true;
data_to_plot = num2cell(abs(analysis_resdyn.result.final_model{1}.eigVal),2);
plotpars.line_style = repmat({'-'},[numel(data_to_plot),1]);
plotpars.marker_style = repmat('o',[length(data_to_plot) 1]);
plotpars.markerfacecol = repmat([1 1 1],[length(data_to_plot) 1]);
[ah] = plotXvsTime(ah,data_to_plot,time_abs,plotpars);
[ah] = formatXvsTimeplot(ah,[],plotpars.do_plot_stability);
yyaxis right
taus_to_show = [-0.1 -0.15 -0.2 -0.35 -0.5 -1 Inf 0.5];
yticks_tau = exp((bin_size./1000)./taus_to_show);
set(ah,'YColor',[0 0 0],'ylim',plotpars.ylim,'ytick',yticks_tau)
y_ticks_tau = arrayfun( @(x) sprintf('%.2f',x), taus_to_show,'uniformoutput', false);
set(ah,'yticklabel',y_ticks_tau);
set(ah,'Parent',gcf);

ah = subplot(1,3,2);set(ah,'plotbox',[1 1 1]);hold(ah)
clear plotpars
clear data_to_plot
plotpars.col_lines = col_lines(round(linspace(1, size(col_lines,1) - 10, size(analysis_resdyn.result.final_model{1}.A,1))),:);
plotpars.do_plot_markers = true;
plotpars.markersize = 4;
plotpars.linewidth = 2;
plotpars.ylim = [-0.5 0.5];
plotpars.ytick = [-0.5 : 0.25 : 0.5];
plotpars.yticklabel = [-0.5 : 0.25 : 0.5];
plotpars.yLabel = 'Angle (eigen-value) (a.u)';
plotpars.do_plot_stability = false;
plotpars.xLabel = 'time (s)';
plotpars.do_plot_stability = false;
data_to_plot = num2cell(analysis_resdyn.result.final_model{1}.eigAngle,2);
plotpars.line_style = repmat({'-'},[numel(data_to_plot),1]);
plotpars.marker_style = repmat('o',[length(data_to_plot) 1]);
plotpars.markerfacecol = repmat([1 1 1],[length(data_to_plot) 1]);
[ah] = plotXvsTime(ah,data_to_plot,time_abs,plotpars);
[ah] = formatXvsTimeplot(ah,[],plotpars.do_plot_stability); 


yyaxis right
set(ah,'YColor',[0 0 0],'ylim',plotpars.ylim,'ytick',plotpars.ytick)
f_to_show = [-2 -1 -0.5 -0.25 0 0.25 0.5 1 2];
yticks_f = 2*pi*(bin_size/1000)*f_to_show;
set(ah,'YColor',[0 0 0],'ylim',plotpars.ylim,'ytick',yticks_f)
y_ticks_f = arrayfun( @(x) sprintf('%.2f',x), f_to_show ,'uniformoutput', false);
set(ah,'yticklabel',y_ticks_f);
set(ah,'Parent',gcf);
 
ah = subplot(1,3,3);set(ah,'plotbox',[1 1 1]);hold(ah)

clear plotpars
clear data_to_plot

plotpars.col_lines = col_lines(round(linspace(1, size(col_lines,1) - 10, size(analysis_resdyn.result.final_model{1}.A,1))),:);
plotpars.do_plot_markers = true;
plotpars.markersize = 4;
plotpars.linewidth = 2;
plotpars.ylim = [0 2.5];
plotpars.ytick = [0 : 0.2 : 2.5];
plotpars.yticklabel = [0 : 0.2 : 2.5];
plotpars.yLabel = '|singular-value | (a.u)';
plotpars.xLabel = 'time (s)';
plotpars.do_plot_stability = true;
data_to_plot = num2cell(abs(analysis_resdyn.result.final_model{1}.sVal),2);
plotpars.line_style = repmat({'-'},[numel(data_to_plot),1]);
plotpars.marker_style = repmat('o',[length(data_to_plot) 1]);
plotpars.markerfacecol = repmat([1 1 1],[length(data_to_plot) 1]);
[ah] = plotXvsTime(ah,data_to_plot,time_abs,plotpars);
[ah] = formatXvsTimeplot(ah,[],plotpars.do_plot_stability);

if(do_save_fig)
    export_fig(fullfile(save_fig_path, sprintf('%s%d%s%.1f%s%.1f','ResDynFits_augmentedlineatt_planeangle=',...
    angle_between_planes,'_lambda_decay=',lambda_decay_1,'_',lambda_decay_2)),'-painters','-transparent','-pdf')
end


%% Computing alingment of cond-avg spaces with dynamics

[data_agg] =  compute_ovelap_resdyn_avgtaskplanes({analysis_condavg}, {analysis_resdyn.result});

%% Saving results to disk

if(do_orth_obs_mat)
    result.data_identifier = sprintf('%s%s%d%s%d','augmentedlineattr','_orthC=',double(do_orth_obs_mat),...
        '_anglebwplanes=',angle_between_planes);
else
    result.data_identifier = sprintf('%s%s%d','augmentedlineattr','_orthC=',double(do_orth_obs_mat));
end

result.resdynfit = data_agg;
result.cond_avgs = {analysis_condavg};
result.resdynfit_pars = opts_analysis;
result.resdynfit_cv = analysis_resdyn;
result.cond_avg_pars.jPC = jPC_pars;
result.cond_avg_pars.choice_time = choice_time_pars;
result.sim_opts = sim_opts;

file_name = strcat(result.data_identifier,'_resdynresults','.mat');
save(fullfile(save_file_path,file_name),'result','-v7.3')

