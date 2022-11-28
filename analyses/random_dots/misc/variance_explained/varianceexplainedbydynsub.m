%{This script computes the varaince explained in the condition-average 
% trajectories by the estimated dynamics subspace. Specifically, this is
% done in two different ways:
% 1) We compute how much of the variance in the 'aligned' choice-averaged
%    neural trajectories unfolding within the individual task planes 
%    (choice, time , jPCxx) areattributable to the dynamcis subspace. This 
%    essentially boils down to computing an overlap between the dynamics
%    subspace and each of the 4  distinct task planes. We also compute a 
%    null dsitribution, by randomly sampling 2D planes in the aligned space
% 2) We compute how much variance is explained in the 'non-aligned',
%    choice-averaged neural trajectories of individual experiemnts by the 
%    dynamics subspace (See Extended Data Fig 7) 
%}
clearvars -except DIRS
clc
rseed = rng('default');
close all
cd('/Users/Aniruddh/Work_PhD/Residual Dynamics/residual_dynamics')

%% Loading all relevant data

animal = 'Tex';
bin_size = 45;
data_file_path = './data/processedneuraldata/';

%Loading resdyn results
file_name = sprintf('%s%s',animal,'_ndim=8_lag=3_alpha=200_50_allconfigsresdynresults.mat');
load(fullfile(data_file_path, file_name));

% Loading aligned single trials
file_name = sprintf('%s%s%d%s',animal,'_aligned_reppasdotsTask_binsize=',bin_size,'ms.mat');
load(fullfile(data_file_path, file_name));

% Loading non-aligned single trials
file_name = sprintf('%s%s%d%s',animal,'_nonaligned_reppasdotsTask_binsize=',bin_size,'ms.mat');
D = load(fullfile(data_file_path, file_name));
data_nonaligned = D.data_aligned;
clear D

% properties
save_fig_path = '';
save_fig_format= '-pdf';
do_save_fig = false;

n_configs = length(data_aligned.aligned_single_trials);
alignments = result.cond_avg_pars.choice_time.alignments;


%% Computing overlap between dynamics subspace and the aligned task planes (along with null)

cols = [0 0 1; 1 0 0];
varexp_8D_2D = cell(1, length(alignments));
varexp_8D_2D_shuf = cell(1, length(alignments));
n_random_dir_sets = 5000;
c_count = 1;
varexp_prc95_flag = [];
varexp_prc99_flag = [];
varexp_orig_vals = [];
n_dim = result.resdynfit_pars.dim{1};

for iexp = 1: n_configs
    figure; set(gcf,'Position',[85  132 1023 573]);
    for iev = 1: length(alignments)     
        
        direction_labels = result.cond_avgs{iexp}.(alignments{iev}).direction_labels;
        n_planes = length(direction_labels); 
        
        for iplanes = 1: length(direction_labels)
            
            ah = subplot(length(alignments), length(direction_labels), (iev - 1)*length(direction_labels) + iplanes);
            hold(ah);set(ah,'plotbox',[1 1 1]);
            
            idx_dir = ismember(result.cond_avgs{iexp}.(alignments{iev}).direction_labels,...
                direction_labels{iplanes});
 
            
            Xmasked_sub = bsxfun(@minus,result.cond_avgs{iexp}.(alignments{iev}).Xavg,...
                mean(result.cond_avgs{iexp}.(alignments{iev}).Xavg(:,:),2));

            U_plane =  result.cond_avgs{iexp}.(alignments{iev}).directions{idx_dir};
            U_dyn = result.resdynfit_final_raw{iexp}.U_dyn(:,1:n_dim);
            
            v1_orig =  trace(cov((U_plane * (U_plane' * Xmasked_sub(:,:)))'));
            v2_orig = trace(cov((U_dyn * (U_dyn' * U_plane) * (U_plane' * Xmasked_sub(:,:)))'));
            varexp_8D_2D{iev}(iplanes) = v2_orig/ v1_orig;
            
            for i_r = 1: n_random_dir_sets
               
                U_rand = myqr((randn(size(U_plane,1) , 2)));
                v1_shuf =  trace(cov((U_rand * (U_plane' * Xmasked_sub(:,:)))'));
                v2_shuf = trace(cov((U_dyn * (U_dyn' * U_rand) * (U_plane' * Xmasked_sub(:,:)))'));        
                varexp_8D_2D_shuf{iev}(iplanes, i_r) = v2_shuf/ v1_shuf;
                
            end
            
            histogram(ah, varexp_8D_2D_shuf{iev}(iplanes,:),'Normalization','pdf');
            prcts = prctile(varexp_8D_2D_shuf{iev}(iplanes,:),[95,99]);
            ylims = get(ah,'ylim');
            plot([varexp_8D_2D{iev}(iplanes) varexp_8D_2D{iev}(iplanes)],ylims,'--k');
            
            varexp_orig_vals(c_count) = varexp_8D_2D{iev}(iplanes);
            if(varexp_8D_2D{iev}(iplanes) >= prcts(1) & varexp_8D_2D{iev}(iplanes) < prcts(2))
                varexp_prc95_flag(c_count) = true;
                varexp_prc99_flag(c_count) = false;
            elseif(varexp_8D_2D{iev}(iplanes) >= prcts(2))
                varexp_prc95_flag(c_count) = true;
                varexp_prc99_flag(c_count) = true;
            elseif(varexp_8D_2D{iev}(iplanes) < prcts(1))
                varexp_prc95_flag(c_count) = false;
                varexp_prc99_flag(c_count) = false;
            end
            c_count = c_count + 1;
            
            set(ah,'xlim',[0 1],'xtick',[0:0.2:1],'xticklabel',[0:0.2:1]);
            if(iev == 2)
                xlabel(ah, 'var exp fraction (v2/v1)')
            end
            
            if(iplanes == 1)
                ylabel(ah, 'probability density (a.u)')
            end
            
            if(iev == 1)
                title(ah,sprintf('%s%s%d%s%s%s%s',animal,'-Cfg',iexp,': '  , ...
                    direction_labels{idx_dir},'-',alignments{iev}));
            end
           
            if(iev == 1 & iplanes == 1)
                legend(ah,{'shuffled','original'});
            end
        end
        
    end
    
    if(do_save_fig)
        figname = sprintf('%s%s%s%d',animal,'_fractiontaskplanevarexpexplainedbydynSub','_config=',iexp);
        export_fig(fullfile(save_fig_path, figname),'-painters','-transparent','-pdf')
    end
end
num_planes_above_95prcs = sum(varexp_prc95_flag)./length(varexp_prc95_flag);
num_planes_above_99prcs = sum(varexp_prc99_flag)./length(varexp_prc99_flag);
%% Computing variance explained by dynamics subspace in non-aligned choice-averaged trajectories of individual experiments.
varexp_dynsub_neuralspace = NaN(1, length(data_nonaligned.dspInfo));
c_count = 1;
for icfg = 1:n_configs
    n_exps = length(data_nonaligned.aligned_single_trials{icfg});
    for iexp = 1:n_exps
        
        [varexp_dynsub_neuralspace(c_count)] = ...
            compute_neuralvariance_explained_by_dynsub(data_nonaligned.aligned_single_trials{icfg}{iexp},...
            result.resdynfit_final_raw{icfg},data_aligned.aligned_single_trials{icfg}.U_align{iexp});
        
        c_count = c_count + 1;
    end
    
end


figure;ah = gca;hold(ah);set(ah,'plotbox',[1 1 1]);
n_bins = 25;
hist_bin_edges = linspace(0,1,25);
plot_col = [0 0 1; 1 0 0];

histogram(unique(varexp_dynsub_neuralspace(:)), hist_bin_edges,'EdgeColor',...
    plot_col(1,:),'FaceColor',plot_col(1,:),'FaceAlpha',0.2)
median_dynsub_varexp = median(unique(varexp_dynsub_neuralspace(:)));
plot([median_dynsub_varexp median_dynsub_varexp],[0 65],'--k');

xlabel('variance explained in neural space (cond-avgs) by dyn subspace')
ylabel('number of experiments');