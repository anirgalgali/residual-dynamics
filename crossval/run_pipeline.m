function [analysis] = run_pipeline(X, time, time_rel, time_labels, trial_labels, opts, rseed)
%{ Function that runs the dynamics fitting pipeline (scluding the session alignment).
% This script has the option of running all the cross-validations and also
% fitting the final model if you want to do everything end-to-end, although
% this is not advised, as you may want to inspect the intermediate outputs
% of the pipeline to determine the optimal values of the hyperparameters.
% Input
% - X (3d tensor)- of size n_dim x n_times x n_trials - Typically the
% residual data matrix.
% - time (array) - of size n_times x 1, indicates the absolute times
%   for all time-bins in the data (typically concatenated over all temporal
%   alignments)
% - time_rel (array) - of size n_times x 1, indicates the 'relative' times
%   in each temporal alignment.
% - time_labels (array) - of size n_times x 1, each element is a
%   numerical assignment of the time bin to a specific task epoch.
% - trial_labels (array) - of size n_trials x 1, each element is a
%   numerical assignment of the trial to a given task condition (for e.g
%   choice)
% - opts(struct) - containing all the options inolved in fitting and
%   cross-validation
% - rseed - matlab random seed
% Output
% - analysis (struct) - containing the output of all the steps of the
%   pipeline (including cv stats, and final fits)
%
% Author - Aniruddh Galgali (First written Jul 2018, Heavily Modified May 2020)
%}
%% Compute variability metrics
unique_trial_labels = unique(trial_labels);
n_conds = length(unique_trial_labels);

unique_time_labels = unique(time_labels);
n_align = length(unique_time_labels);

for icond = 1: n_conds
    for ialign = 1:n_align
        [analysis.var_metrics.l0(icond,ialign),analysis.var_metrics.l1(icond,ialign),...
            analysis.var_metrics.l1_l0(icond,ialign),analysis.var_metrics.p_var(icond,ialign)]...
            = compute_prop_ratio(X(:,time_labels == unique_time_labels(ialign),trial_labels == unique_trial_labels(icond)));
    end
end
analysis.fit_pars = opts;
analysis.fit_pars.rseed = rseed;

%% Hankel rank cross-validation

hankel_order = opts.hankel_order;

if(size(X,1) > 1)
    
    hankel_ranks_cv = opts.hankel_cv.hankel_ranks;
    n_cv_h = opts.hankel_cv.n_cv;
    criterion = opts.hankel_cv.criterion;
    metric_type = opts.hankel_cv.type;
    
    
    cv_hankel = cell(n_conds,1);
    dims_out = cell(n_conds,1);
    opt_hankel_rank = cell(n_conds,1);
    
    
    for icond = 1:n_conds
        
        [cv_hankel{icond}, dims_out{icond}] = cross_validate_hankelrank(X(:,:,trial_labels == ...
            unique_trial_labels(icond)),hankel_order, hankel_ranks_cv, time, time_labels, n_cv_h, rseed);
        
        assert(any(ismember(fieldnames(dims_out{icond}),criterion)),'invalid criterion')
        assert(any(ismember(fieldnames(dims_out{icond}.(criterion)),metric_type)),'invalid metric')
        
        opt_hankel_rank{icond} = dims_out{icond}.(criterion).(metric_type);
        
        
    end
else
   opt_hankel_rank = repmat({1},[n_conds,1]);
end

analysis.cv_hankel = cv_hankel;
analysis.cv_hankel_opt = dims_out;
analysis.cv_hankel_hyp = hankel_ranks_cv;
analysis.opt_hankel_rank = opt_hankel_rank;


%% Run Lag DIM CV

[analysis.cv_lagdim, analysis.cv_lagdim_opt, analysis.cv_lagdim_hyp] = ...
    hyperparamsearch_laganddim(X, hankel_order, opt_hankel_rank,...
    time, trial_labels, time_labels, rseed, opts.lag_cv);

%%

opt_lag = NaN(n_conds ,1);
opt_dim = NaN(n_conds ,1);

for icond = 1: n_conds
    
    switch opts.final_fit.lag_cv_metric
        
        case 'mse'
            
            opt_lag(icond) = analysis.cv_lagdim_opt{icond}.mse_opt.lags(1);
            opt_dim(icond) = analysis.cv_lagdim_opt{icond}.mse_opt.dims(1);
            
        case 'pvar'
            
            opt_lag(icond) = analysis.cv_lagdim_opt{icond}.pvar_opt.lags(1);
            opt_dim(icond) = analysis.cv_lagdim_opt{icond}.pvar_opt.dims(1);
    end
    
    
    
end

if(opts.final_fit.use_common_lagdim_across_conditions)
    num_pars = opt_lag.*(opt_dim.^2);
    [~,idx] = min(num_pars);
    opt_lag = opt_lag(idx).*ones(n_conds,1);
    opt_dim = opt_dim(idx).*ones(n_conds,1);
end

%% run smoothness cv

if(opts.doSmoothCV)
     
    assert(isfield(opts,'smooth_cv'),'field containing smoothing hyperparameters is missing');
    [U_dyn_full] = compute_dynamics_subspace(X,trial_labels,time,time_labels,hankel_order,opt_hankel_rank);
    [cv_smoothness, cv_smooth_opt] = hyperparamsearch_alpha(X,U_dyn_full,opt_lag,opt_dim,trial_labels,time_labels, opts.smooth_cv,rseed);
    
    analysis.cv_smoothness = cv_smoothness;
    analysis.cv_smooth_opt = cv_smooth_opt;
    analysis.cv_smooth_hyp = opts.smooth_cv.alpha_all;
    analysis.cv_smooth_pars.lag = opt_lag;
    analysis.cv_smooth_pars.dim = opt_dim;

end

%%
if(opts.doFinalFit)
    
    if(opts.doSmoothCV)
            
        alpha_all_mse = cell2mat(cellfun(@(x) x.T_mse.alpha(end), analysis.cv_smooth_opt,'uni',false));
        alpha_all_pvar = cell2mat(cellfun(@(x) x.T_pvar.alpha(end), analysis.cv_smooth_opt,'uni',false));
        
        for icond = 1:n_conds
            switch opts.final_fit.alpha_cv_metric
                
                case 'pvar'
                    opt_final_fit.alpha{unique_trial_labels(icond)} = alpha_all_pvar(:,unique_trial_labels(icond));
                    
                case 'mse'
                    opt_final_fit.alpha{unique_trial_labels(icond)} = alpha_all_mse(:,unique_trial_labels(icond));
                    
            end
        end

    
    else
        
        assert(isfield(opts.final_fit,'smoothalpha'),'need to provide a valid regulairzation parameter')
        assert(length(opts.final_fit.smoothalpha) == n_align,'number of alphas should match number of alignments');
        opt_final_fit.alpha = repmat(opts.final_fit.smoothalpha,[n_conds,1]);

        
    end
    opt_final_fit.lag = num2cell(opt_lag,2);
    opt_final_fit.dim = num2cell(opt_dim,2);
    opt_final_fit.hankel_order = opts.hankel_order;
    opt_final_fit.hankel_rank = analysis.opt_hankel_rank;
    opt_final_fit.do_cross_validate_noise_model = opts.final_fit.do_cross_validate_noise_model;
    
    [analysis.result] = fit_residual_dynamics(X, time, time_rel, time_labels, trial_labels, opt_final_fit);
    

end


end


