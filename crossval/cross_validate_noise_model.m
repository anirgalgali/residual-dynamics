function [err_iso_avg,err_aniso_avg,err_tv_avg,err_tv_iso_avg] = ...
    cross_validate_noise_model(X,U_dyn, opt_lag, opt_alpha, time_labels)

n_folds = 5;
cv_inds = crossvalind('KFold',size(X,3), n_folds);

unique_time_labels = unique(time_labels);
nalign = length(unique_time_labels);

% err_iso_avg = zeros(n_folds,nalign);
% err_aniso_avg = zeros(n_folds,nalign);
% err_tv_avg = zeros(n_folds,nalign);
% err_tv_iso_avg = zeros(n_folds,nalign);

for ii = 1:n_folds
   
    test_inds = cv_inds == ii;
    train_inds = cv_inds ~= ii;
    
    [data_fit] = applydynamicsSubspaceProjection(X(:,:,train_inds), U_dyn);
    [A, time_idxs_dyn, X_fs] = estimateDynamicsTwoStageLS(data_fit, opt_lag, opt_alpha, time_labels);
    time_labels_res = time_labels(time_idxs_dyn);
   
    for ialign = 1:nalign
        
        Xev = X(:,time_labels == unique_time_labels(ialign),:);
        
        [R, Q0, Q_iso_fixed, Q_aniso_fixed, Q_tv, Q_tv_iso] = ...
        computeLatentCovarianceParametersMLE(A(:,:,time_labels_res == unique_time_labels(ialign)), U_dyn, Xev(:,opt_lag + 1 : end, train_inds), X_fs{unique_time_labels(ialign)});
    
    
        X_test = Xev(:,opt_lag + 1:end, test_inds);
        X_test = applydynamicsSubspaceProjection(X_test, U_dyn);
        X_test_cov_iso = NaN(size(X_test,1),size(X_test,1),size(X_test,2) - 1);
        X_test_cov_aniso = NaN(size(X_test,1),size(X_test,1),size(X_test,2) - 1);
        X_test_cov_tv = NaN(size(X_test,1),size(X_test,1),size(X_test,2) - 1);
        X_test_cov_tv_iso = NaN(size(X_test,1),size(X_test,1),size(X_test,2) - 1);
        
        err_iso = zeros(size(X_test,2)-1,1);
        err_aniso = zeros(size(X_test,2)-1,1);
        err_tv = zeros(size(X_test,2)-1,1);
        err_tv_iso = zeros(size(X_test,2)-1,1);
        
        for tt = 1: size(X_test,2)-1
            
            X_test_cov_iso(:,:,tt) = A(:,:,tt)*cov(squeeze(X_test(:,tt,:))')*A(:,:,tt)' + Q_iso_fixed(:,:,tt);
            X_test_cov_aniso(:,:,tt) = A(:,:,tt)*cov(squeeze(X_test(:,tt,:))')*A(:,:,tt)' + Q_aniso_fixed(:,:,tt);
            X_test_cov_tv(:,:,tt) = A(:,:,tt)*cov(squeeze(X_test(:,tt,:))')*A(:,:,tt)' + Q_tv(:,:,tt);
            X_test_cov_tv_iso(:,:,tt) = A(:,:,tt)*cov(squeeze(X_test(:,tt,:))')*A(:,:,tt)' + Q_tv_iso(:,:,tt);
            
            err_iso(tt) = norm(cov(squeeze(X_test(:,tt+1,:))') - X_test_cov_iso(:,:,tt) ,'fro').^2;
            err_aniso(tt) = norm(cov(squeeze(X_test(:,tt+1,:))') - X_test_cov_aniso(:,:,tt) ,'fro').^2;
            err_tv(tt) = norm(cov(squeeze(X_test(:,tt+1,:))') - X_test_cov_tv(:,:,tt) ,'fro').^2;
            err_tv_iso(tt) = norm(cov(squeeze(X_test(:,tt+1,:))') - X_test_cov_tv_iso(:,:,tt) ,'fro').^2;
            
        end
        
        err_iso_avg(ii,:,ialign) = err_iso;
        err_aniso_avg(ii,:,ialign) = err_aniso;
        err_tv_avg(ii,:,ialign) = err_tv;
        err_tv_iso_avg(ii,:,ialign) = err_tv_iso;
    end
    
end


end