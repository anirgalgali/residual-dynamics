%{This script plots the various plots connected with the session alignment
% procedure (see Extended Data Figure 6c-g) for a single animal.
%
% Author: Aniruddh Galgali (Modified: May 2021)
%}

clearvars -except DIRS
clc
rseed = rng('default');
close all

%% Choose monkey and dataset to analuze
animal = 'Tex';
bin_size = 45;
do_save_fig = false;
save_fig_path = '' ;% provide a valid path for saving the figures;
%% Load aligned data

cd('/Users/Aniruddh/Work_PhD/Residual Dynamics/residual_dynamics')
data_path = './data/processedneuraldata/';
file_name = sprintf('%s%s%d%s',animal,'_aligned_reppasdotsTask_binsize=',bin_size,'ms.mat');
load(fullfile(data_path,file_name)); 
%% Extracting and plotting aligned condition-averaged responses for a single configuration

% Choosing the configuration to plot

n_configs = length(data_aligned.aligned_single_trials);
config_idx_to_plot = 3; % (config_idx_to_plot < n_configs)
[n_modes, n_times, n_conds, n_exps] = size(data_aligned.alignment_info.projs{config_idx_to_plot});
time_labels = data_aligned.aligned_single_trials{config_idx_to_plot}.time_iev;
time_abs_all = data_aligned.aligned_single_trials{config_idx_to_plot}.time;
time_rel_all = data_aligned.aligned_single_trials{config_idx_to_plot}.time_rel;
unique_time_labels = unique(data_aligned.aligned_single_trials{config_idx_to_plot}.time_iev);
n_align = length(unique_time_labels);

% extracting aligned projections
choice_averaged_data = num2cell(permute(data_aligned.alignment_info.projs{config_idx_to_plot},[1 2 4 3]),[1 2 3]);

%Plotting average aligned projections as line plots (Supplementary Figure 5b)

figure; set(gcf,'Position',[151 613 1018 429]);
if(config_idx_to_plot == 1)
    col_lines = {'-r';'-b'}; % choice indices for config-1 are flipped!
else
    col_lines = {'-b';'-r'};
end
for imode = 1:n_modes
   
   ah = subplot(n_modes/5,n_modes/4,imode);hold(ah);set(ah,'plotbox',[1 1 1]);
   
   for iev = 1: length(unique_time_labels)
       time_idxs = time_labels == unique_time_labels(iev);
       for icond = 1:n_conds
           
           aligned_projection_mean = mean(choice_averaged_data{icond}(imode,time_idxs,:),3);
           aligned_projection_sem = (1/sqrt(size(choice_averaged_data{icond},3))).*...
               std(choice_averaged_data{icond}(imode,time_idxs,:),0,3);
           
           shadedErrorBar(time_abs_all(time_idxs),aligned_projection_mean,...
               2.*aligned_projection_sem,col_lines{icond})

           
       end
 
   end
   
   set(gca,'ylim',[-3 3],'ytick',[-3 0 3],'yticklabel',[-3 0 3]);
   if(imode == 1)
       xlabel(gca,'time(s)','FontSize',8)
       ylabel(gca,'projection (a.u)','FontSize',8)
   end
   
    
end
suptitle(sprintf('%s%s%d',animal,', Config-',config_idx_to_plot));
if(do_save_fig)
    save_fig_format = '-pdf';
    fig_name = sprintf('%s%s',id_list_selected{iexp},'_alignedProjs_Mean&Sem_allsessions');
    export_fig(fullfile(save_fig_path,fig_name),'-painters','-transparent',save_fig_format);
end

% Plotting aligned projs for all experiemnts in a single configuration as a
% single image plot.

n_plots = n_modes *n_conds*n_align;
n_rows = (n_modes/5)*n_conds;
n_cols = (n_modes/4)*n_align;

% This is a hacky way of indicating ther order in which the different
% modes (column1), choices(cilumn2) and alignment (column3) are plotted.
plot_indices = NaN(n_plots,3);
plot_indices(:,2) = [ones(n_cols,1);2*ones(n_cols,1);ones(n_cols,1);2*ones(n_cols,1);...
    ones(n_cols,1);2*ones(n_cols,1);ones(n_cols,1);2*ones(n_cols,1)];
plot_indices(:,3) = repmat([1;2],[n_plots/2 1]);
plot_indices(:,1) = [1;1;2;2;3;3;4;4;5;5;1;1;2;2;3;3;4;4;5;5;...
                     6;6;7;7;8;8;9;9;10;10;6;6;7;7;8;8;9;9;10;10;...
                     11;11;12;12;13;13;14;14;15;15;11;11;12;12;13;13;14;14;15;15;...
                     16;16;17;17;18;18;19;19;20;20;16;16;17;17;18;18;19;19;20;20];


figure; set(gcf,'Position',[151 613 1018 429]);

for iplot = 1:n_plots
    
    ah = subplot((n_modes/5)*n_conds,(n_modes/4)*n_align,iplot);
    hold(ah);set(ah,'plotbox',[1 1 1]);hold all;
    time_idxs = (time_labels == plot_indices(iplot,3));
    hh = imagesc(squeeze(choice_averaged_data{plot_indices(iplot,2)}(plot_indices(iplot,1),time_idxs,:))');
    
    if(plot_indices(iplot,1) == 1)
        cc = 4;
    elseif(plot_indices(iplot,1) > 1 & plot_indices(iplot,1) <= 4)
        cc = 2;
    elseif(plot_indices(iplot,1) >4 & plot_indices(iplot,1) <= 7)
        cc = 1;
    else
        cc = 0.5;
    end
    caxis([-cc cc]);
    
    axis tight;colormap('jet');
    if(plot_indices(iplot,2) == 1 & plot_indices(iplot,3) == 1)
        
        if(plot_indices(iplot,1) == 1 || plot_indices(iplot,1) == 2 || plot_indices(iplot,1) == 5 || plot_indices(iplot,1) == 8)
            cbh = colorbar;
            cbh.Ticks = [-cc, 0 , cc];
            cbh.TickLabels = num2cell([round(-cc,1), 0 , round(cc,1)]);
            cbh.Location = 'north';
        end
    end
    set(hh,'Parent',ah);
    
    time_rel = time_rel_all(time_idxs);
    
    if(plot_indices(iplot,3) == 1)
        
        [~,idx1] = min(abs(time_rel - 0));
        [~,idx2] = min(abs(time_rel - 0.8));
        plot([idx1 idx1],[1 size(choice_averaged_data{plot_indices(iplot,2)},3)],'-k','linewidth',0.2)
        plot([idx2 idx2],[1 size(choice_averaged_data{plot_indices(iplot,2)},3)],'-k','linewidth',0.2)
        
    elseif(plot_indices(iplot,3) == 2)
        
        [~,idx1] = min(abs(time_rel - 0));
        plot([idx1 idx1],[1 size(choice_averaged_data{plot_indices(iplot,2)},3)],'-k','linewidth',0.2)
    end
    
    
    if(plot_indices(iplot,3) == 2)
        set(gca,'ytick',[],'yticklabel',[]);
    else
        set(gca,'ytick',[10:10:40],'yticklabel',[10:10:40]);
    end
    set(gca,'xtick',[],'xticklabel',[]);
    
    if(iplot == 1)
        xlabel(gca,'time (s)','FontSize',8)
        ylabel(gca,'exp #','FontSize',8);
    end
    
end
suptitle(sprintf('%s%s%d',animal,', Config-',config_idx_to_plot));
if(do_save_fig)
    fig_name = sprintf('%s%s',id_list_selected{iexp},'_alignedProjs_allsessions_implotMaxColrange');
    export_fig(fullfile(save_fig_path,fig_name),'-painters','-transparent','-pdf');
end


%% Plotting aligned correlation matrices for all configurations (Supplementary Figure 5e)

figure;
for iconfig = 1:n_configs 
    
    corrMat_median = data_aligned.alignment_info.stats{iconfig}.align_sim;
    ah = subplot(1, n_configs , iconfig);hold(ah);set(ah,'plotbox',[1 1 1]);
    hh = imagesc(corrMat_median);axis tight;colormap(flipud(cbrewer('div','RdYlBu',64)));
    caxis([0 1])
    set(hh,'Parent',ah);
    set(ah,'YDir','normal');
    title(ah,sprintf('%s%s%d',animal,', Config-',iconfig))
    xlabel(ah,'Mode #')
    ylabel(ah,'Mode #')
end

if(do_save_fig)
    save_fig_format = '-pdf';
    fig_name = sprintf('%s%s','SessionAlignCorrMat','_allConfigs');
    export_fig(fullfile(save_fig_path,fig_name),'-painters','-transparent',save_fig_format);
end

%%  Plotting session alignment variance explained

config_markers = {'^';'d';'o';'s'};
n_modes_var_plot = 100;
n_mode_var_markers = [1 5 10:10:100]; 
n_modes = 20;
figure;set(gcf,'Position',[680  452 1048 646]);
ah1 = subplot(1,3,1); hold(ah1); set(ah1,'plotbox',[1 1 1]);
ah2 = subplot(1,3,2); hold(ah2); set(ah2,'plotbox',[1 1 1]);
ah3 = subplot(1,3,3); hold(ah3); set(ah3,'plotbox',[1 1 1]);

for iconfig = 1: n_configs
    
    var_explained_config = data_aligned.alignment_info.stats{iconfig}.varExp_ses;
    var_explained_config = cell2mat(var_explained_config');
    var_explained_config = var_explained_config(:,1:n_modes_var_plot);
    
    n_exps = size(var_explained_config,1);
    num_units_per_exp = cell2mat(cellfun(@(x) size(x,1), data_aligned.aligned_single_trials{iconfig}.U_align,'uni',false));
    num_trials_per_exp = NaN(n_exps,1);
    for iexp = 1:n_exps
        num_trials_per_exp(iexp) = sum(data_aligned.aligned_single_trials{iconfig}.task_variable.session_id == iexp);
    end
    
    plot(ah1,1 : n_modes_var_plot, mean(var_explained_config,1),'-','Color',[0 0 0]);
    errorbar(ah1,n_mode_var_markers , mean(var_explained_config(:,n_mode_var_markers),1),...
        (1/sqrt(n_exps))./std(var_explained_config(:,n_mode_var_markers),0,1),...
        config_markers{iconfig},'Color',[0 0 0],'MarkerSize',6,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 1 1]);
    
    plot(ah2,num_units_per_exp,var_explained_config(:,n_modes),...
        config_markers{iconfig},'MarkerSize',3,'MarkerEdgeColor',[0 0 0],...
        'MarkerFaceColor',[0 0 0]);

    plot(ah3,log10(num_trials_per_exp), var_explained_config(:,n_modes),...
        config_markers{iconfig},'MarkerSize',3,'MarkerEdgeColor',[0 0 0],...
        'MarkerFaceColor',[0 0 0]);

end
set(ah1,'xlim',[0 100],'xtick',n_mode_var_markers,'xticklabel',n_mode_var_markers);
set(ah1,'ylim',[0 100],'ytick',[0:20:100],'yticklabel',[0:20:100]);
plot(ah1,[n_modes n_modes],[0 100],'--','Color',[0.5 0.5 0.5]);
xlabel(ah1,'# of aligned modes')
ylabel(ah1,'variance explained (%)')
xlabel(ah2,'# of units in individual exps')
ylabel(ah2,'variance explained (%)')
xlabel(ah3,'# of trials in individual exps')
ylabel(ah3,'variance explained (%)')
suptitle(sprintf('%s%s',animal,': Session alignment'))

%%
% Loading the cross-validation results, which containg information about
% U_dyn
file_name = sprintf('%s%s%d%s',animal,'_aligned_hankelagdimCV_reppasdotsTask_binsize=',bin_size,'ms.mat');
R = load(fullfile(data_path,file_name));
result_lagdimCV = R.result;

n_dynsub_dim = 8;
marker_cols = cbrewer('div','Spectral',n_dynsub_dim);
config_markers = {'^';'d';'o';'s'};
marker_size = 3;
figure;ah = gca;hold(ah);set(ah,'plotbox',[1 1 1]);
do_xaxis_noise = true;

for iconfig = 1 : n_configs
    for idim = 1:n_dynsub_dim
        if(~do_xaxis_noise)
            plot(1:n_modes, abs(result_lagdimCV{iconfig}.U_dyn(:,idim)),config_markers{iconfig},...
                'markeredgecolor', marker_cols(idim,:),'markerfacecolor',marker_cols(idim,:),'markersize',2);
        else
            plot([1:n_modes] + 0.1*randn(1,n_modes), abs(result_lagdimCV{iconfig}.U_dyn(:,idim)),config_markers{iconfig},...
                'markeredgecolor', marker_cols(idim,:),'markerfacecolor',marker_cols(idim,:),'markersize',marker_size);
            
        end
        set(ah,'ylim',[0 1],'ytick',[0:0.2:1],'yticklabel',[0:0.2:1]','xlim',[0 21]);
        
    end
end
xlabel('aligned mode')
ylabel('|projection|')
title(sprintf('%s%s',animal,' - alignment of dyn subspace with aligned subspace'))
if(do_save_fig)
    save_fig_format = '-pdf';
    fig_name = sprintf('%s%s%d','DynamicsSubspaceLoadings','_allConfigs_ndynsubdim=',n_dynsub_dim);
    export_fig(fullfile(save_fig_path,fig_name),'-painters','-transparent',save_fig_format);
end