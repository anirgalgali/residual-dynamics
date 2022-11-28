function [outA,varargout] = bootstrapDynamics(X,iv_lag,lambda,time_labels,n_boot,n_time_dyn,type, A_stats)

[n_dim,n_time,~] = size(X);
X = permute(X,[3 1 2]);
X = X(:,:);
boot_opts = statset;
boot_opts.UseParallel = true;

switch type
    
    case 'direct'
        
        boot_fn = @(x) bootstrapTVLDS(x,iv_lag,lambda,time_labels,n_dim,n_time);
        [bstat,bootsam] = bootstrp(n_boot,boot_fn,X,'Options',boot_opts);
        A_boot = NaN(n_dim,n_dim,n_time_dyn,n_boot);
        outA = cell(n_boot,1);
        
        for i_boot = 1:n_boot
            A_boot(:,:,:,i_boot) = reshape(bstat(i_boot,:)',[n_dim,n_dim,n_time_dyn]);
            outA{i_boot}.A = A_boot(:,:,:,i_boot);
        end
        
    case 'ci'
        
        boot_fn = @(x) bootstrapTVLDS_bootci(x,iv_lag,lambda,time_labels,n_dim,n_time, A_stats);
        [ci_out, bstat] = bootci(n_boot,{boot_fn,X},'alpha',0.05,'Options',boot_opts);
        
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