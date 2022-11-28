function [cv_out,dims_out,varargout] = cross_validate_hankelrank(X,order,ranks,time,time_labels,n_cv,rand_state)
%{ Function that runs the cross-validation of the hankel ranks.
% Input
% - X (3d tensor)- of size n_dim x n_times x n_trials - Typically the
% residual data matrix.
% - order (scalar) - the hankel order (q)
% - ranks (array) - specifying all the ranks (hyperparams) that you want to
%   test as part of the cross-validation.
% - time (array) - size n_times x 1 , which indicates the time-bins
%    associated with the data. Typically, this is an array constructed by
%    concatanating times across all temporal alignments. 
% - time_labels (array) - size n_times x 1 , each entry (numerical) indicates
%   the label indicating the assignment of the time bin to a task epoch.
% - n_cv(scalar) - number of repeats of hold-oout cross validation
% - rand_state - MATLAB random seed.
% Output
% - cv_out (struct) - containing the cross-validated fit errors and the
%   associated statsistics
% - dims_out (struct) - containing the optimal hankel ranks for different
%   fit error criteria.
%
% Author - Aniruddh Galgali (First written Jul 2018)
%}
rng(rand_state);
n_trials = size(X,3);
n_ranks = length(ranks);

train_recons_error = cell(n_cv,1);
test_recons_error = cell(n_cv,1);

for i_cv = 1:n_cv

    cv_inds = crossvalind('HoldOut',n_trials);
    data_train = X(:,:,cv_inds == 0);
    data_test = X(:,:,cv_inds == 1);
    
    Xtrain = data_train - repmat(mean(data_train,3),[1 1 size(data_train,3)]);
    Xtest = data_test - repmat(mean(data_test,3),[1 1 size(data_test,3)]);

    [H_train] = computeTimeVaryingHankelMatrix(Xtrain, time, order, time_labels);
    [H_test,time_h] = computeTimeVaryingHankelMatrix(Xtest, time, order, time_labels);
    
    train_error = NaN(size(H_train,3),n_ranks);
    test_error = NaN(size(H_train,3),n_ranks);
    
    for jj = 1: n_ranks
        
        [H_lr] = hankelLowRankApprox(H_train,ranks(jj));
        
        for tt = 1:size(H_lr,3)
            
            train_error(tt,jj) = norm(H_train(:,:,tt) - H_lr(:,:,tt),'fro');
            test_error(tt,jj) = norm(H_test(:,:,tt) - H_lr(:,:,tt),'fro');
            
            
        end
        
    end
    
    train_recons_error{i_cv} = train_error;
    test_recons_error{i_cv} = test_error;
    
end


idxs_time_h = ismember(time,time_h);
time_labels_h = time_labels(idxs_time_h);
[cv_out, dims_out] = computeCVstats(train_recons_error,test_recons_error,time_labels_h, ranks);
varargout{1} = time_h;
varargout{2} = time_labels_h;

end

function [cv_out,dims_out] = computeCVstats(train_error,test_error,time_labels,ranks)

train_error = cat(3,train_error{:});
test_error = cat(3,test_error{:});
unique_time_labels = unique(time_labels);
train_error_timeAveraged_combined = squeeze(mean(train_error,1));
test_error_timeAveraged_combined = squeeze(mean(test_error,1));

train_error_timeAveraged_align = cell(length(unique_time_labels),1);
test_error_timeAveraged_align = cell(length(unique_time_labels),1);
for ii = 1 : length(unique_time_labels)
    
    train_error_timeAveraged_align{ii} = squeeze(mean(train_error(time_labels == unique_time_labels(ii),:,:),1));
    test_error_timeAveraged_align{ii} = squeeze(mean(test_error(time_labels == unique_time_labels(ii),:,:),1));
end

cv_out.train_error = train_error;
cv_out.test_error = test_error;
cv_out.ranks = ranks;

cv_out.stats.train_mean.combined = mean(train_error_timeAveraged_combined,2);
cv_out.stats.train_std.combined = std(train_error_timeAveraged_combined,0,2);
cv_out.stats.test_mean.combined = mean(test_error_timeAveraged_combined,2);
cv_out.stats.test_std.combined = std(test_error_timeAveraged_combined,0,2);

cv_out.stats.train_mean.align = cellfun(@(x) mean(x,2),train_error_timeAveraged_align,'uni',false);
cv_out.stats.train_std.align = cellfun(@(x) std(x,0,2),train_error_timeAveraged_align,'uni',false);
cv_out.stats.test_mean.align = cellfun(@(x) mean(x,2),test_error_timeAveraged_align,'uni',false);
cv_out.stats.test_std.align = cellfun(@(x) std(x,0,2),test_error_timeAveraged_align,'uni',false);

cv_out.stats.train_mean.time_varying = mean(train_error,3);
cv_out.stats.train_std.time_varying  = std(train_error,0,3);
cv_out.stats.test_mean.time_varying  = mean(test_error,3);
cv_out.stats.test_std.time_varying  = std(test_error,0,3);

[~,idx_c] = min(cv_out.stats.test_mean.combined);
[~,idx_a] = cellfun(@(x) min(x),cv_out.stats.test_mean.align,'uni',false);
[~,idx_t] = min(cv_out.stats.test_mean.time_varying,[],2);


dims_out.min.combined = ranks(idx_c);

for jj = 1:length(cv_out.stats.test_mean.align)
    dims_out.min.align(jj) = ranks(idx_a{jj});
end

for jj = 1:size(cv_out.stats.test_mean.time_varying,1)
    dims_out.min.time_varying(jj) = ranks(idx_t(jj));
end

uB_combined = cv_out.stats.test_mean.combined(idx_c) + ...
    (1/sqrt(size(train_error,3))).*cv_out.stats.test_std.combined(idx_c);

valid_ranks = ranks(cv_out.stats.test_mean.combined <= uB_combined);
dims_out.sem.combined = min(valid_ranks);

uB_combined = cv_out.stats.test_mean.combined(idx_c) + ...
    cv_out.stats.test_std.combined(idx_c);

valid_ranks = ranks(cv_out.stats.test_mean.combined <= uB_combined);
dims_out.std.combined = min(valid_ranks);

for jj = 1:length(cv_out.stats.test_mean.align)
    
    uB_align = cv_out.stats.test_mean.align{jj}(idx_a{jj}) + ...
        (1/sqrt(size(train_error,3))).*cv_out.stats.test_std.align{jj}(idx_a{jj});
    valid_ranks = ranks(cv_out.stats.test_mean.align{jj} <= uB_align);
    dims_out.sem.align(jj) = min(valid_ranks);
    
    
    uB_align = cv_out.stats.test_mean.align{jj}(idx_a{jj}) + ...
        cv_out.stats.test_std.align{jj}(idx_a{jj});
    valid_ranks = ranks(cv_out.stats.test_mean.align{jj} <= uB_align);
    dims_out.std.align(jj) = min(valid_ranks);
    
    
end

for jj = 1:size(cv_out.stats.test_mean.time_varying,1)
    uB_tV = cv_out.stats.test_mean.time_varying(jj,idx_t(jj)) + ...
        (1/sqrt(size(train_error,3))).*cv_out.stats.test_std.time_varying(jj,idx_t(jj));
    valid_ranks = ranks(cv_out.stats.test_mean.time_varying(jj,:) <= uB_tV);
    dims_out.sem.time_varying(jj) = min(valid_ranks);
    
    
    uB_tV = cv_out.stats.test_mean.time_varying(jj,idx_t(jj)) + ...
        cv_out.stats.test_std.time_varying(jj,idx_t(jj));
    valid_ranks = ranks(cv_out.stats.test_mean.time_varying(jj,:) <= uB_tV);
    dims_out.std.time_varying(jj) = min(valid_ranks);
    
end

end

