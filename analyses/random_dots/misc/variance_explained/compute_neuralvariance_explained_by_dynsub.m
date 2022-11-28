function [varexp_dyn_all, varargout] = compute_neuralvariance_explained_by_dynsub(data_st, result_final, Ualign)
% This function computes the variance explained by the dynamics subspace in
% the condition-avered trajectories in the full neural space.  
% Input
% - data_st (struct) - contains the non-aligned single trials and associated 
%   variables for a particular experiment
% - result_final (struct) - contains the fits of the residual dynamics for
%   a specific task configuration. 
% - Ualign (matrix) - of size n_units x n_align - which is the orthogonal
%   projection matrix the relates the high-dim neural space to the aligned
%   space.
% Output
% - varexp_dyn_all (scalar) - the variance explained

%  When arguments are passed to this function, make sure that 
%  data_st and result final belong to the same task configuration for the
%  results to make sense.
%  Author: Aniruddh Galgali (April 2022).

n_rand_dirs = 5000;

trial_labels = data_st.task_index.targ_dir;
unique_time_labels = unique(data_st.time_iev);
unique_cond_labels = unique(trial_labels);

n_dim = size(result_final.final_model{1}.A,1);
Xcond = [];
for iev = 1: length(unique_time_labels)
    
    X = data_st.response(:,data_st.time_iev == unique_time_labels(iev),:);
    Xcond_ = NaN(size(X,1), size(X,2), length(unique_cond_labels));
    
    for icond = 1 : length(unique_cond_labels)
        
        Xcond_(:,:,unique_cond_labels(icond)) = mean(X(:,:,trial_labels == unique_cond_labels(icond)),3);
        
    end
    Xcond = cat(2,Xcond,Xcond_);
end
    
Xcond_msub = bsxfun(@minus, Xcond, mean(Xcond,2));
origVar = trace(cov(Xcond_msub(:,:)'));
proj_dynsub = Ualign * result_final.U_dyn(:,1:n_dim);
projVar_dynSub = trace(cov((proj_dynsub* (proj_dynsub'*Xcond_msub(:,:)))'));
varexp_dyn_all = projVar_dynSub./origVar;

[varexp_rand_dynsub_constrained] = compute_varexp_randdirs_dynsubconstrained(Xcond_msub,Ualign,...
     result_final.U_dyn(:,1:n_dim), n_rand_dirs);
  

varargout{1} = varexp_rand_dynsub_constrained;

end


function [varexp] = compute_varexp_randdirs_dynsubconstrained(X,U_align,U_dyn,n_rand_dirs)
    
    origVar = trace(cov(X(:,:)'));
    rand_wts =  randn(size(U_dyn,2), n_rand_dirs);
    rand_dirs = normc(U_align*(U_dyn*rand_wts));
    varexp = NaN(1,n_rand_dirs);

    for idir = 1: n_rand_dirs
        
        projVar = trace(cov((rand_dirs(:,idir)* (rand_dirs(:,idir)'*X(:,:)))'));
        varexp(idir) = projVar/origVar;
    end

end

