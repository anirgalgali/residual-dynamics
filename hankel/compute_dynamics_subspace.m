function [U,varargout] = compute_dynamics_subspace(X,trial_labels,time,time_labels,h_order,h_rank)

%{ 

Computes a time-invariant dynamics subspace U from an intermediate
sequence of time-varying hankel matrices H(t) that are constructed using
the data in the 3D data tensor X. 

Inputs : X -            3D data tensor of size n_dim x n_time x n_trials
         trial_labels - array of size n_trials x1 indicating the condition
                        label associated with each trial. The computation of
                        the intermediate  time-varying hankel matrices is 
                        done separately for each condition.
         time -         array of size n_time x 1 indicating the time stamps
                        associated with data in X. The time axis here is a 
                        single contiguous one across all alignments
         time_labels -  array of size n_time x 1, each point labeling the 
                        event alignment associated with a given time bin in
                        time
         h_order -      scalar specifying the order used to compute the
                        hankel matrices
         h_rank -       cell array of length n_conds x 1, where n_conds
                        indicates the number of conditions (i.e number of
                        unique elements in trial_labels). Each element in
                        h_rank can either be :
                        
                        (i)   a scalar - specifies a single rank used for 
                              computing the time-varying observation 
                              matrices C_t from the intermediate sequence of 
                              hankel matrices H_t. OR
                        (ii)  an array of size n_align x 1, where n_align is
                              the number of event alignments (i.e number of
                              unique elements in time_labels). In this 
                              scenario, each element of this array specifies 
                              the rank used for computing the time-varying 
                              observation matrices C_t for a particular event 
                              alignment.   OR
                        (iii) an array of size n_times_hankel x 1, where
                              n_times_hankel is the sequence length of the
                              time-varying hankel matrices. In this
                              scenario, each time-varying observation
                              matrix C_t is computed using a corresponding
                              rank that is stored in this array.
                    

Outputs : U - matrix of size n_dim x n_dim, where the columns are ordered in
             decreasing order of importance
          OPTIONAL (in-order):
            1) time_var_obsmats - A 

        Author : Aniruddh Galgali, Dec 2017
                 Major modification: May 2020

%}
n_trials = size(X,3);
assert(n_trials == length(trial_labels),'number of trials do not match');
unique_cond_labels = unique(trial_labels);
n_conds = length(unique_cond_labels);


C_t_all = cell(1,n_conds);
cond_labels_Ct_all = cell(1,n_conds);
time_Ct_all = cell(1,n_conds);
time_labels_Ct_all = cell(1,n_conds);

for icond = 1:n_conds

    X_fit = X(:,:,trial_labels == unique_cond_labels(icond));
    X_fit = X_fit - repmat(mean(X_fit,3),[1 1 size(X_fit,3)]);    % ensure that mean across trials is zero
    [C_t, time_hankel, time_label_hankel] = ...
    compute_timevarying_observation_matrices(X_fit, time, time_labels, h_order, h_rank{icond});

    C_t_all{icond} = C_t;
    cond_labels_Ct_all{icond} = unique_cond_labels(icond).*ones(1,length(C_t));
    time_Ct_all{icond} = time_hankel;
    time_labels_Ct_all{icond} = time_label_hankel;
    
end
C_t_all = cat(2,C_t_all{:});
cond_labels_Ct_all = cell2mat(cond_labels_Ct_all);
time_Ct_all = cell2mat(time_Ct_all);
time_labels_Ct_all = cell2mat(time_labels_Ct_all);

obsmats_t.C = C_t_all;
obsmats_t.cond_labels = cond_labels_Ct_all;
obsmats_t.time = time_Ct_all;
obsmats_t.time_labels = time_labels_Ct_all;

[U,S,~] = svd(cell2mat(C_t_all),'econ');   

[obsmats_t.overlap_with_dynsub] = computeoverlap_Ct_and_dynamicssubspace(C_t_all, U);

varargout{1} = obsmats_t;
varargout{2} = S;

end


function [C_t, varargout] = compute_timevarying_observation_matrices(X, time, time_labels, hankel_order, hankel_rank)

%{ 
This function is a wrapper that computes the time-varying observation matrices 
from the time-varying hankel matrices. 
%}
n_dim = size(X,1);
unique_time_labels = unique(time_labels);
num_align = length(unique_time_labels);

[H_t,time_hankel,time_labels_hankel] = computeTimeVaryingHankelMatrix(X, time, hankel_order, time_labels);
n_times_hankel = size(H_t,3);
    
if(numel(hankel_rank) == 1)
    
    opt_rank = repmat(hankel_rank,[1 n_times_hankel ]);
    
elseif(numel(hankel_rank) == num_align)
    
    opt_rank = [];
    unique_time_labels_hankel = unique(time_labels_hankel);
    assert(length(unique_time_labels_hankel) == num_align,'wrong size');
    
    for ialign = 1:num_align
        
        n_times_hankel_align = sum(time_labels_hankel == unique_time_labels_hankel(ialign));
        opt_rank = [opt_rank ...
            repmat(hankel_rank(unique_time_labels_hankel(ialign)), [1 n_times_hankel_align])];
        
    end
    
elseif(numel(hankel_rank) == n_times_hankel)
    
    opt_rank = hankel_rank;
    
else
    
    error('array containing ranks is of wrong size');
end
    
[O_t] = computeObservabilityMatrixfromHankel(H_t, opt_rank);
C_t = cellfun(@(x) x(1 : n_dim , :), O_t, 'uni', false);

varargout{1} = time_hankel;
varargout{2} = time_labels_hankel;

end