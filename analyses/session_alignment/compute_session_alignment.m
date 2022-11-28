function [D_align,varargout] = compute_session_alignment(data,condition_vars,pars)
%{ Performs the session alignment
%  Input
%  - data (cell array of structs) - of size n_experiments x 1. Each element 
%    of the cell array is a struct that contains the single trial data
%    belonging to a single experiment
%  - condition_vars (cell array) - specifies the task variables based on 
%    which the trial-averaged trajectories are computed for session
%    alignment.
%  - pars (struct) - struct containing the high-level parameters that
%    determine the session alignment
%  Output
%  - D_align (cell array of structs) - similar to data, but contains the
%    aligned single-trial responses.
%  - varargout
%    - D_cond_avg_proj_align (cell array of structs) - simislar to D_align,
%    but contains only contains the trial-averaged projections of the 
%    responses used for plotting
%    -align_stats (struct) - structure containing the statistics and
%     information about the alignment.
%
%  Author: Aniruddh Galgali (May 2018, Modified: Jul 2020)
%
%}
% check input data types
assert(iscell(data),'variable: data - must be a cell array');
assert(iscell(condition_vars),'variable: task_conditions - must be a cell array');

% Pre-processing datasets - binning/session exclusions/...
num_sessions = length(data);

% Sorting configuration based datasets into task conditions for alignment
% computations
   
row_count = 0;
row_idxs_all = cell(1,num_sessions);
D_align = cell(1, num_sessions);
varExp_ses = cell(1, num_sessions);
D_avg = [];

for i_session =  1: num_sessions
    [D_cond] = sort_trials_by_condition(data{i_session},condition_vars);
    num_units = size(D_cond{1}.response,1);
    row_idxs = row_count + 1 : row_count + num_units;
    row_idxs_all{i_session} = row_idxs;
    D_avg(row_idxs,:,:) = cell2mat(cellfun(@(x) mean(x.response,3),D_cond,'uni',false)');
    row_count = row_count + num_units;
end
    
[coeff,~,~,~,explained,mu] = pca(D_avg');
varExp_pc = cumsum(explained');
[U_orth] = orthogonalize_sub_matrices(coeff,row_idxs_all);

for i_session = 1:num_sessions
    
    D_align{i_session} = data{i_session};
    [~,num_times,num_trials] = size(data{i_session}.response);
    meanSub_response_full = (data{i_session}.response(:,:) - mu(row_idxs_all{i_session})');
    aligned_response = U_orth{i_session}(:,1:pars.alignment.nPCs_align)' * meanSub_response_full;      
    D_align{i_session}.response = reshape(aligned_response,[pars.alignment.nPCs_align num_times num_trials]);
    D_align{i_session}.U = U_orth{i_session}(:,1:pars.alignment.nPCs_align);
    D_align{i_session}.X_mu = mu(row_idxs_all{i_session})';
    
    for idim = 1:size(U_orth{i_session},2)
        varExp_ses{i_session}(idim) = compute_session_variance(D_avg(row_idxs_all{i_session},:) - ...
            mu(row_idxs_all{i_session})', U_orth{i_session}(:,1:idim)); 
    end

end
    

D_cond_avg_proj_align = cell(1,num_sessions);
% align_sim_5 = cell(num_configs,1);
% align_sim_95 = cell(num_configs,1);
for i_session =  1: num_sessions     

    [D_cond_align,task_conds_proj] = sort_trials_by_condition(D_align{i_session},pars.align_proj.condition_vars);        
    D_avg_align = cellfun(@(x) mean(x.response,3),D_cond_align,'uni',false);
    D_avg_align = cat(3,D_avg_align{:});
    D_cond_avg_proj_align{i_session} = D_avg_align;


end 
D_cond_avg_proj_align = cat(ndims(D_avg_align) + 1, D_cond_avg_proj_align{:});  
[align_sim_med] = compute_alignment_similarity(D_cond_avg_proj_align);


align_stats.varExp_pc = varExp_pc;
align_stats.varExp_ses = varExp_ses;
align_stats.align_sim = align_sim_med; 
varargout{1} = D_cond_avg_proj_align;
varargout{2} = task_conds_proj;
varargout{3} = align_stats;
end


function [U_orth] = orthogonalize_sub_matrices(U,unit_indices)

num_sessions = length(unit_indices);
U_orth = cell(num_sessions,1);
for i_seg = 1: num_sessions
   
    U_seg = U(unit_indices{i_seg},:);
    [Q,~] = myqr(U_seg(:,1:min(size(U_seg))));
    U_orth{i_seg} = Q;
end

end

function [var_exp] = compute_session_variance(D,U)
var_exp = percvar(D, U*U'*D);
end

function [M_med] = compute_alignment_similarity(D)

n_modes = size(D,1);
M_med = NaN(n_modes,n_modes);

for i1 = 1:n_modes
    for i2 = 1:n_modes
        X = permute(squeeze(D(i1,:,:,:)),[3 1 2]);
        Y = permute(squeeze(D(i2,:,:,:)),[3 1 2]);
        M = corr(X(:,:)',Y(:,:)');
        mask = triu(true(size(M)),1);
        M_med(i1,i2) = median(M(mask));

    end
end

end
