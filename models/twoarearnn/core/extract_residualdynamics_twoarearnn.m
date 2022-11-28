function [result_fit] = extract_residualdynamics_twoarearnn(residuals, time, time_rel,...
                            time_labels,trial_labels, analysis, fit_opts, varargin)
%{ Fits the final residual dynamics to simulated residuals from a single network 
%  configurations
%  
%  Input
%  - residuals - (p x T x K) - 3d tensor containing residuals at each time and trial 
%  - time - (1 x T) - array containing the absolute times corresponding to the
%                     residuals
%  - time_rel - (1 x T) - array containing the 'relative' times (within each epoch)
%                      corresponding to the residuals
%  - time_labels - (1 x T) - array containing the time labels indicating the epoch
%                   for each time
%  - trial_labels - (K x 1) - array containing the condition labels
%                   corresponding to each trial
%  - analysis - struct containing the results of the cross-validaiton of
%               the optiaml hyperparams
%  - fit_opts - additional parameters for the final fits
%  - varargin  - (i) specifies a single state-space direction along which you may
%                   want to compute residual dynamics
%
%  Output
%  - result_fit - struct containing the final fits
%
%  Author: Aniruddh Galgali (Oct 2020)
%}    
minArgs = 7;
maxArgs = 8;
narginchk(minArgs,maxArgs)

if(nargin == 8)
   direction = varargin{1};
else
   direction = [];
end

n_conds = length(analysis.cv_lagdim_opt);

opts = struct();
for icond = 1:n_conds
    if(isfield(fit_opts,'dim') & isempty(direction))
        
        % If the parameters is specified as a field in the fit_opts
        % structure, then the setting of the parameter obtained through grid
        % search is over-ridden (only if you aren't a specifying a specific
        % direction along which to compute the 
        
        opts.dim{icond} = fit_opts.dim;
    elseif(~isfield(fit_opts,'dim') & isempty(direction))
        assert(isfield(fit_opts, 'max_dim'),'you need to provide the ceiling dim')
        opts.dim{icond} = min(analysis.cv_lagdim_opt{icond}.mse_opt.dims(1),fit_opts.max_dim);
    else
        % If a direction is specified, then the residual dynamics is just
        % 1D.
        opts.dim{icond} = 1;
    end
    
    if(isfield(fit_opts,'lag'))
        opts.lag{icond} = fit_opts.lag;
    else
        opts.lag{icond} = analysis.cv_lagdim_opt{icond}.mse_opt.lags(1);
    end
    
    if(isfield(fit_opts,'alpha'))
        opts.alpha{icond} = fit_opts.alpha;
    else
        opts.alpha{icond} = analysis.cv_smooth_opt{icond}.T_mse.alpha(1); % Largest value of alpha
    end
end


if(isempty(direction))
    
   opts.do_boot = false;
   if(isfield(fit_opts,'bootstrp'))
       opts.boot_opts = fit_opts.bootstrp;
   end
   % Hankel order and ranks are only needed to compute the dynamics
   % subspace, which is only necessary if a direction is not specified.
   opts.hankel_order = analysis.fit_pars.hankel_order;
   opts.hankel_rank = cellfun(@(x) x.(fit_opts.hankel_criterion).(fit_opts.hankel_type),analysis.cv_hankel_opt,'uni',false);
   opts.do_cross_validate_noise_model = false;
   [result_fit] = fit_residual_dynamics(residuals, time, time_rel, time_labels, trial_labels, opts);

elseif(~isempty(direction))
    
   opts.do_boot = true;
   if(isfield(fit_opts,'bootstrp'))
       opts.boot_opts = fit_opts.bootstrp;
   else
       error('need to provide bootstrap options');
   end
   opts.do_cross_validate_noise_model = false;
   residuals = applydynamicsSubspaceProjection(residuals, direction);
   [result_fit] = fit_residual_dynamics(residuals, time, time_rel, time_labels, trial_labels, opts);
   
end

end




