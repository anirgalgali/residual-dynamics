%{
This scripts performs all the statistical tests reported in the paper in
relation to the results presented in Figure 4.

Author - Aniruddh Galgali (Apr 2022)
%}

clearvars -except DIRS
close all
clc
%%

monkey = 'Vito';
bin_size = 45;
load_file_path = './data/processedneuraldata';
data_file_name = sprintf('%s%s',monkey,'_ndim=8_lag=3_alpha=200_50_allconfigsresdynresults.mat');
load(fullfile(load_file_path, data_file_name))
dt =  bin_size/1000;

%% Retaining only those entries of the table that correspond to the largest eigenvalue

max_ev_all = result.resdynfit(result.resdynfit.is_max == 1,:);
unique_time_labels = unique(max_ev_all.time_label);

%% Largest EV < 1 ?

ev_filt_cond = {};
pvals = cell(1,length(unique_time_labels));
n = cell(1,length(unique_time_labels),1);    
for iep = 1: length(unique_time_labels)
    
    idx_ep = max_ev_all.time_label == unique_time_labels(iep);
    all_unique_times_in_epoch = unique(max_ev_all.time_rel(idx_ep));
    
    if(~isempty(ev_filt_cond))
        
        idx_times = all_unique_times_in_epoch >= ev_filt_cond{unique_time_label(iep)}(1) & ...
            all_unique_times_in_epoch <= ev_filt_cond{unique_time_label(iep)}(2);
        
        times_to_analyze = all_unique_times_in_epoch(idx_times);
        
    else
        
        times_to_analyze = all_unique_times_in_epoch;
   
    end
    for tt = 1: length(times_to_analyze)

            idx_tt = max_ev_all.time_rel == times_to_analyze(tt);
            max_ev_tt = max_ev_all.ev(idx_tt);
            [h, p] = ttest(max_ev_tt, 1,'Tail','left');
            
            pvals{iep}(tt) = p;
            n{iep}(tt) = length(max_ev_tt);
        
    end
end
tail_prob = [0.005 0.01]; %{0.01 or 0.005}
all_pvals = cell2mat(pvals);
all_n = cell2mat(n);
fprintf('Fraction of points with p < %.3f (n = %d) : %d/%d\n',tail_prob(1),unique(all_n),sum(all_pvals<tail_prob(1)),length(all_pvals));
fprintf('Fraction of points with p < %.3f (n = %d) : %d/%d\n',tail_prob(2),unique(all_n),sum(all_pvals<tail_prob(2)),length(all_pvals));

%% Median time-constants and their CIs

ev_filt_cond = {[0 0.8],[-0.5 0.3]}; % use this for time windows of interest
% ev_filt_cond = {}; % use this if you want the median across all times
taus = cell(1,length(unique_time_labels));

epoch_names = {'stim';'sacc'};
for iep = 1: length(unique_time_labels)
    
    idx_ep = max_ev_all.time_label == unique_time_labels(iep);
    all_unique_times_in_epoch = unique(max_ev_all.time_rel(idx_ep));
    
    if(~isempty(ev_filt_cond))
        
        idx_times = all_unique_times_in_epoch >= ev_filt_cond{unique_time_labels(iep)}(1) & ...
            all_unique_times_in_epoch <= ev_filt_cond{unique_time_labels(iep)}(2);
        
        times_to_analyze = all_unique_times_in_epoch(idx_times);
        
    else
        
        times_to_analyze = all_unique_times_in_epoch;
   
    end
    for tt = 1: length(times_to_analyze)

            idx_tt = max_ev_all.time_rel == times_to_analyze(tt);
            max_ev_tt = max_ev_all.ev(idx_tt);
            tau_max_ev_tt = dt./log(max_ev_tt);
            taus{iep}(:,tt) = tau_max_ev_tt;
            
        
    end
    n = length(taus{iep}(:));
    median_tau = median(taus{iep}(:));
    cis = prctile(taus{iep}(:),[2.5 97.5]);
    fprintf('%s : Median tau = %.4f, 95th percentile CIs : %.4f , %.4f (n = %d)\n',epoch_names{iep},median_tau,cis(1),cis(2),n);

end
n = length(vec(cell2mat(taus)));
median_tau = median(vec(cell2mat(taus)));
cis = prctile(vec(cell2mat(taus)),[2.5 97.5]);
fprintf('All times: Median tau = %.4f, 95th percentile CIs : %.4f , %.4f (n = %d)\n',median_tau,cis(1),cis(2),n);

%% Comparing maximum eigenvalues in delay versus peri-saccade times

probe_time_delay = -0.275;
probe_time_perisacc = -0.005;
max_ev_delay= max_ev_all(max_ev_all.time_rel == probe_time_delay & max_ev_all.time_label == 2,:);
max_ev_perisacc= max_ev_all(max_ev_all.time_rel == probe_time_perisacc & max_ev_all.time_label == 2,:);
[h, p] = ttest2(max_ev_delay.max_ev, max_ev_perisacc.max_ev,'Tail','right');
fprintf('H0:e.v equal at %.3f s and %.3f s relative to sacc onset,  p = %.6f (n = %d)\n',probe_time_delay,probe_time_perisacc,p,size(max_ev_delay,1));

%% Comparing max_ev veruss max_sv

ev_filt_cond = {};
pvals = cell(1,length(unique_time_labels));
n = cell(1,length(unique_time_labels),1);    
for iep = 1: length(unique_time_labels)
    
    idx_ep = max_ev_all.time_label == unique_time_labels(iep);
    all_unique_times_in_epoch = unique(max_ev_all.time_rel(idx_ep));
    
    if(~isempty(ev_filt_cond))
        
        idx_times = all_unique_times_in_epoch >= ev_filt_cond{unique_time_label(iep)}(1) & ...
            all_unique_times_in_epoch <= ev_filt_cond{unique_time_label(iep)}(2);
        
        times_to_analyze = all_unique_times_in_epoch(idx_times);
        
    else
        
        times_to_analyze = all_unique_times_in_epoch;
   
    end
    for tt = 1: length(times_to_analyze)

            idx_tt = max_ev_all.time_rel == times_to_analyze(tt);
            max_ev_tt = max_ev_all.ev(idx_tt);
            max_sv_tt = max_ev_all.max_sv(idx_tt);
            [h, p] = ttest2(max_sv_tt,max_ev_tt,'Tail','right');
            
            pvals{iep}(tt) = p;
            n{iep}(tt) = length(max_ev_tt);
        
    end
end
tail_prob = [0.05]; %{0.01 or 0.005}
all_pvals = cell2mat(pvals);
all_n = cell2mat(n);
fprintf('Fraction of points with p < %.3f (n = %d) : %d/%d\n',tail_prob(1),unique(all_n),sum(all_pvals<tail_prob(1)),length(all_pvals));
%%