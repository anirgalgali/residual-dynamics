%{
This scripts plots the aligned neural activity for each task configuration 
within 4 distinct "task-planes". These task planes capture variance due to
choice, time and rotational dynamics (jPCxx). See Figure 3
Author - Aniruddh Galgali (Nov 2020)
%}

clearvars -except DIRS
close all
clc
%%

monkey = 'Tex';
bin_size = 45;
load_file_path = './data/processedneuraldata';
data_file_name = sprintf('%s%s',monkey,'_ndim=8_lag=3_alpha=200_50_allconfigsresdynresults.mat');
load(fullfile(load_file_path, data_file_name))
dt =  bin_size/1000;
do_save_fig = false;
save_fig_path = ''; % enter a valid path;

%%
n_configs = length(result.cond_avgs);
config_idx_to_plot = 3;
alignments = result.cond_avg_pars.choice_time.alignments;
direction_labels = result.cond_avgs{1}.(alignments{1}).direction_labels;
pref_line_col = [0 0 1]; 
anti_line_col = [1 0 0];
axis_pad = 0.2;
pref_idx = result.cond_avg_pars.choice_time.pref_idx;
anti_idx = result.cond_avg_pars.choice_time.anti_idx;

for icfg = 1:n_configs
    
    figure;set(gcf,'Position',[680         295        1131         803]);

    for ialign = 1: length(alignments)
        for idir = 1: length(direction_labels)

            idx_dir = ismember(result.cond_avgs{icfg}.(alignments{ialign}).direction_labels, direction_labels{idir});

            projection_to_plot = result.cond_avgs{icfg}.(alignments{ialign}).projections{idx_dir};
            ah = subplot(length(alignments),length(direction_labels),(ialign-1)*length(direction_labels) + idir);
            hold(ah);set(ah,'plotbox',[1 1 1]);

            plot(squeeze(projection_to_plot(1,:,pref_idx(icfg))),squeeze(projection_to_plot(2,:,pref_idx(icfg))),...
                '-','Color',pref_line_col);
            plot(squeeze(projection_to_plot(1,:,anti_idx(icfg))),squeeze(projection_to_plot(2,:,anti_idx(icfg))),...
                '-','Color',anti_line_col);

            plot(squeeze(projection_to_plot(1,:,pref_idx(icfg))),squeeze(projection_to_plot(2,:,pref_idx(icfg))),...
                'o','markeredgecolor',pref_line_col,'markerfacecolor',[1 1 1],'markersize',4);
            plot(squeeze(projection_to_plot(1,:,anti_idx(icfg))),squeeze(projection_to_plot(2,:,anti_idx(icfg))),...
                'o','markeredgecolor',anti_line_col,'markerfacecolor',[1 1 1],'markersize',4);

            [~,ev_idx] = min(abs(result.cond_avgs{icfg}.(alignments{ialign}).time - 0));

            plot(squeeze(projection_to_plot(1,ev_idx,pref_idx(icfg))),squeeze(projection_to_plot(2,ev_idx,pref_idx(icfg))),...
                'o','markeredgecolor',pref_line_col,'markerfacecolor',pref_line_col,'markersize',4);
            plot(squeeze(projection_to_plot(1,ev_idx,anti_idx(icfg))),squeeze(projection_to_plot(2,ev_idx,anti_idx(icfg))),...
                'o','markeredgecolor',anti_line_col,'markerfacecolor',anti_line_col,'markersize',4);

            axis tight
            ylims = get(gca,'ylim');
            xlims = get(gca,'xlim');
            yrange = ylims(2) - ylims(1);
            xrange = xlims(2) - xlims(1);
            max_range = max(xrange,yrange);
            max_range = max_range + axis_pad;
            set(ah,'ylim',[(ylims(2) + ylims(1))/2 - max_range/2 (ylims(2) + ylims(1))/2 + max_range/2]);
            set(ah,'xlim',[(xlims(2) + xlims(1))/2 - max_range/2 (xlims(2) + xlims(1))/2 + max_range/2]);
            title(ah,sprintf('%s%s%s',direction_labels{idir}, ' ,',alignments{ialign}))
            
            
            if(idir == 1 & ialign == 1)
               xlabel(ah,'dim-1 (a.u)')
               ylabel(ah,'dim-2 (a.u)')
            end
            
        end
        
    end
    
    if(do_save_fig)
        fig_name = sprintf('%s%s%d%s',animal,'Config=',icfg,'_TimeChoice&jPCProjections');
        export_fig(fullfile(save_fig_path,fig_name),'-painters','-transparent','-pdf');
    end

end