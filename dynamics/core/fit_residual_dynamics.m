function [result] = fit_residual_dynamics(X, time, time_rel, time_labels, trial_labels, opts)

% Input Arguments:
% X - residual data matrix of size p x T x K (p : aligned dim, T : number
% of time points, K : numner of trials)

% time - time vector specifying the binned time value corresponding to each
% data point in X

% time_labels - a categorical array that specifies the label corresponding
% to each time point in X. The label specifies the alignment of the data in
% X and determines how the data in X are split before computing the
% residual dynamics

% trial_labels - a categorical array that specifies the label corresponding
% to each trial. This labels determines how the trials in X are grouped in
% order to compute the residual dynamics
% 
% opts (struct) - containing the parameters that determine the fit (such as
%                 the hyperparams, and other options)
% 
% Output
%
% result(struct) - containging the time-varying residual dynamics and all
%                  its properties
%
% Author: Aniruddh Galgali (Apr 2019, modified in May 2020)

assert(length(trial_labels) == size(X,3),'data matrix is of incorrect size');
assert(length(time_labels) == size(X,2),'data matrix is of incorrect size');

unique_trial_labels = unique(trial_labels);
n_conds = length(unique_trial_labels);
unique_time_labels = unique(time_labels);
n_align = length(unique_time_labels);

if(isfield(opts,'hankel_order') | isfield(opts,'hankel_rank'))
    
    assert(length(opts.hankel_rank) == n_conds,...
        'number of ranks to input should equal number of residual conditions');
    [result.U_dyn, result.C_t] = compute_dynamics_subspace(X, trial_labels, time, time_labels,...     
              opts.hankel_order, opts.hankel_rank);
end
         
result.fit_pars = opts;      
for icond = 1: n_conds
    

    
    assert(length(opts.dim{icond}) == 1,...
        'both alignments should use the same latent dimensionality');
    assert(length(opts.lag{icond}) == 1,...
        'both alignments should use the same lag');
    assert(length(opts.alpha{icond}) == n_align);
    
    if(isfield(result, 'U_dyn'))
        [data_fit] = applydynamicsSubspaceProjection(X(:,:,trial_labels == unique_trial_labels(icond)),...
            result.U_dyn(:,1:opts.dim{icond}));
    elseif(~isfield(result, 'U_dyn') & size(X,1) == 1)
        data_fit = X(:,:,trial_labels == unique_trial_labels(icond));
    else
        error('need to provide a dynamics subspace');
    end
    
    [A,time_idxs_dyn,X_fs] = estimateDynamicsTwoStageLS(data_fit,...
        opts.lag{icond}, opts.alpha{icond}, time_labels);
    
    result.final_model{icond}.A = A;
    
    if(isfield(opts, 'do_save_fs_trajs'))
        if(opts.do_save_fs_trajs)
            result.final_model{icond}.X_fs = X_fs;
        end
    end
    
    result.final_model{icond}.time = time(time_idxs_dyn);
    result.final_model{icond}.time_iev = time_labels(time_idxs_dyn);
    result.final_model{icond}.time_rel = time_rel(time_idxs_dyn);
    
    if(isfield(opts,'order_across_alignment'))
        [outA] = computeDerivedDynamicsQuantities(A,'abs','right',time_labels(time_idxs_dyn), time(2) - time(1), opts.order_across_alignment);
    else
        [outA] = computeDerivedDynamicsQuantities(A,'abs','right',time_labels(time_idxs_dyn), time(2) - time(1));

    end
    outAfnames = fieldnames(outA);
    for ff = 1: length(outAfnames)
        
        result.final_model{icond}.(outAfnames{ff}) = outA.(outAfnames{ff});
        
    end
    
    if(opts.do_cross_validate_noise_model)
        
        
        [result.final_model{icond}.cv_noise.err_iso_avg,result.final_model{icond}.cv_noise.err_aniso_avg,...
            result.final_model{icond}.cv_noise.err_tv_avg,result.final_model{icond}.cv_noise.err_tv_iso_avg] = ...
            cross_validate_noise_model(X(:,:,trial_labels == unique_trial_labels(icond)),...
            result.U_dyn(:,1:opts.dim{icond}), opts.lag{icond}, opts.alpha{icond}, time_labels);
        

    end
    
    if(isfield(opts, 'do_save_fs_trajs'))
        for ialign = 1:n_align
            Xev = X(:,time_labels == unique_time_labels(ialign),:);
            [result.final_model{icond}.noiseCovs{unique_time_labels(ialign)}.R, result.final_model{icond}.noiseCovs{unique_time_labels(ialign)}.Q0,...
                result.final_model{icond}.noiseCovs{unique_time_labels(ialign)}.Q_iso_fixed, result.final_model{icond}.noiseCovs{unique_time_labels(ialign)}.Q_aniso_fixed,...
                result.final_model{icond}.noiseCovs{unique_time_labels(ialign)}.Q_aniso_tv, result.final_model{icond}.noiseCovs{unique_time_labels(ialign)}.Q_iso_tv] = ...
                computeLatentCovarianceParametersMLE(result.final_model{icond}.A(:,:,result.final_model{icond}.time_iev == unique_time_labels(ialign)), result.U_dyn(:,1:opts.dim{icond}), ...
                Xev(:,opts.lag{icond}+1:end,trial_labels == unique_trial_labels(icond)),result.final_model{icond}.X_fs{unique_time_labels(ialign)});
        end
    end
    
    if(isfield(opts,'do_boot'))
       if(opts.do_boot)
           [result.final_model{icond}.boot_stats] = bootstrap_dynamics(data_fit, opts.lag{icond}, opts.alpha{icond},...
               time_labels, result.final_model{icond}, outA, opts.boot_opts);
       end
    end
  
end

end

function [outA, varargout] = bootstrap_dynamics(X, iv_lag, lambda, time_labels, true_result, A_stats, b_opts)

[n_dim,n_time,~] = size(X);
X = permute(X,[3 1 2]);
X = X(:,:);
boot_opts = statset;
boot_opts.UseParallel = b_opts.use_parallel;

switch b_opts.type
    
    case 'direct'
        
        boot_fn = @(x) bootstrapTVLDS(x,iv_lag,lambda,time_labels,n_dim,n_time);
        [bstat,bootsam] = bootstrp(b_opts.n_boot,boot_fn,X,'Options',boot_opts);
        A_boot = NaN(n_dim,n_dim,length(true_result.time),b_opts.n_boot);
        outA = cell(b_opts.n_boot,1);
        
        for i_boot = 1:b_opts.n_boot
            A_boot(:,:,:,i_boot) = reshape(bstat(i_boot,:)',[n_dim,n_dim,length(true_result.time)]);
            [outA{i_boot}] = align_bootstrap(A_stats,A_boot(:,:,:,i_boot), true_result.time(2) - true_result.time(1));
            outA{i_boot}.A = A_boot(:,:,:,i_boot);
        end
        
    case 'ci'
        
        boot_fn = @(x) bootstrapTVLDS_bootci(x,iv_lag,lambda,time_labels,n_dim,n_time, A_stats);
        [ci_out, bstat] = bootci(b_opts.n_boot,{boot_fn,X},'alpha',0.05,'Options',boot_opts);
        
        outA.eigVal = ci_out(:,1:n_dim,:);
        outA.eigAngle = ci_out(:,n_dim+1:2*n_dim,:);
        outA.Re_lambda = ci_out(:,2*n_dim+1:3*n_dim,:);
        outA.Im_lambda = ci_out(:,3*n_dim+1:4*n_dim,:);
        outA.sVal = ci_out(:,4*n_dim+1:5*n_dim,:);
        outA.sVal_sq = ci_out(:,5*n_dim+1:6*n_dim,:);
        outA.traceCov = ci_out(:,6*n_dim+1,:);
        outA.depN = ci_out(:,end-1,:);
        outA.max_EV = ci_out(:,end,:);
        bootsam =[];
end
varargout{1} = bstat;
varargout{2} = bootsam;

end


function [A] = bootstrapTVLDS(X,iv_lag,lambda,time_labels,n_dim,n_time)

assert(ismatrix(X),'need data to be reshaped into a 2d matrix instead of tensor');
n_trials = size(X,1);
X = permute(reshape(X,[n_trials n_dim n_time]),[2 3 1]);
[A] = estimateDynamicsTwoStageLS(X,iv_lag,lambda,time_labels);

end

function [out_all] = bootstrapTVLDS_bootci(X,iv_lag,lambda,time_labels,n_dim,n_time, A_stats)

assert(ismatrix(X),'need data to be reshaped into a 2d matrix instead of tensor');
n_trials = size(X,1);
X = permute(reshape(X,[n_trials n_dim n_time]),[2 3 1]);
[A] = estimateDynamicsTwoStageLS(X,iv_lag,lambda,time_labels);

[junk, D_ev_boot] = my_eigenshuffle_pairwise(A_stats.eigSpec, A_stats.eigVec, A, A_stats.ev_sort_type);
[junk1, D_sv_boot,junk2] = my_singularshuffle_pairwise(A_stats.sVal, A_stats.rSVec, A_stats.lSVec, A, A_stats.sv_sort_type);

D_sv_squared = D_sv_boot.^2;
traceCov = sum(D_sv_squared,1);
depN = sqrt(sum(D_sv_squared,1) - sum(abs(D_ev_boot).^2,1))./sqrt(sum(abs(D_ev_boot).^2,1));

out_all = cat(1,abs(D_ev_boot),angle(D_ev_boot),...
    real(D_ev_boot),imag(D_ev_boot), D_sv_boot,...
    D_sv_squared, traceCov, depN, max(abs(D_ev_boot),[],1));


end

function [out_boot] = align_bootstrap(A_ref_stats,A_boot, dt)

[V_ev_boot, D_ev_boot] = my_eigenshuffle_pairwise(A_ref_stats.eigSpec, A_ref_stats.eigVec,A_boot, A_ref_stats.ev_sort_type);
[U_sv_boot, D_sv_boot,V_sv_boot] = my_singularshuffle_pairwise(A_ref_stats.sVal, A_ref_stats.rSVec, A_ref_stats.lSVec, A_boot, A_ref_stats.sv_sort_type);

out_boot.eigVal = abs(D_ev_boot);
out_boot.eigAngle = angle(D_ev_boot);
[out_boot.tau] =  dt./log(dt./log(out_boot.eigVal));
[out_boot.rotFreq] =  out_boot.eigAngle./(2*pi*dt);
out_boot.eigVec = V_ev_boot;
out_boot.eigSpec = D_ev_boot;
out_boot.sVal = D_sv_boot;
out_boot.lSVec = U_sv_boot;
out_boot.rSVec = V_sv_boot;
out_boot.traceCov = sum(D_sv_boot.^2,1);
out_boot.depN = sqrt(out_boot.traceCov  - sum(abs(D_ev_boot).^2,1))./sqrt(sum(abs(D_ev_boot).^2,1));
out_boot.ev_sort_type = A_ref_stats.ev_sort_type;
out_boot.sv_sort_type = A_ref_stats.sv_sort_type;

end
