function [ci_masked] = create_bootci_masks(est,ci,t_idx, time_rel, time_iev)
%{ Helper function to extract the bootstrap CIs and the true estimate of the 
%  eigenvalues at specific time indices
%}
unique_time_labels= unique(time_iev);
assert(length(t_idx) == length(unique_time_labels),'unqueal number of alignments');

ci_masked = NaN(size(ci));

for iEv = 1: length(unique_time_labels)
    
    t_idx_ci_ev = t_idx{iEv};
    time_rel_ev = time_rel(time_iev == unique_time_labels(iEv));
    true_vals = est(:,time_iev == unique_time_labels(iEv));
    ci_ev = ci(:,time_iev == unique_time_labels(iEv),:);
    ci_masked_ev = NaN(size(ci_ev));
    
    for ii = 1: numel(t_idx_ci_ev)
        
        [~,idx_bars] = min(abs(time_rel_ev  - t_idx_ci_ev(ii)));
        [~,idx_dim_sorted] = sort(true_vals(:,idx_bars),'descend');
        idx_dim_ci = [idx_dim_sorted(1) idx_dim_sorted(end)];
        ci_masked_ev(idx_dim_ci, idx_bars, :) = ci_ev(idx_dim_ci, idx_bars, :);
        
    end
    
    ci_masked(:,time_iev == unique_time_labels(iEv),:) = ci_masked_ev;
    
end


end