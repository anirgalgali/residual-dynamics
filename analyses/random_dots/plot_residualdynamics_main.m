%{ This script plots the different components of the main residual dynamics 
%  figure (Figure 4). Specifically, it shows all the ev and sv for a single
%  animal and task configuration, along with loadings of the dynamics
%  subspace onto the aligned space (Fig 4a-d)
%  It also plots the summary of the maximum ev and sv (Fig 4d) for both
%  animals.
%
%  This scipt is a bit of a monstrosity, and can be cleaned to make it look
%  prettier.
%
%  Author - Aniruddh Galgali (Written - May 2021)
%}

clearvars -except DIRS
close all
clc

%% Loading stored results

animals = {'Tex';'Vito'};
bin_size = 45;
load_file_path = './data/analyses/';
result = cell(length(animals),1);
for iani =1:length(animals)
    data_file_name = sprintf('%s%s',animals{iani},'_ndim=8_lag=3_alpha=200_50_allconfigsresdynresults.mat');
    R = load(fullfile(load_file_path, data_file_name));
    result{iani} = R.result;
end

%% Selecting animal and config to show and specifying whether to show bootstrap results

animal = 'Tex';
config_idx_to_plot = 3;
do_plot_bootstrap = true;
cidx = 1; % choice index to plot

if(do_plot_bootstrap)
    
    file_name = sprintf('%s%s%d%s',animal,'_resdynbootstrapped_config',...
            config_idx_to_plot,'_dim=8_lag=3_alpha=200_50.mat');
    result_bootstrap = load(fullfile(load_file_path,file_name));
    result{ismember(animals,animal)}.resdynfit_final_raw{config_idx_to_plot} = result_bootstrap.result_final;
end
%% Setting figure properties

n_rows = 7;
n_cols = 8;
figure;set(gcf,'Position',[243 265 1486 765]);

%% Plotting the loadings of the dynamics subspace (Fig 4a)

ah = subplot(n_rows,n_cols,[1 9]);hold(ah);
dynsub_dim = result{ismember(animals,animal)}.resdynfit_final_raw{config_idx_to_plot}.fit_pars.dim{cidx};
imagesc(1:dynsub_dim,1:size(result{ismember(animals,animal)}.resdynfit_final_raw{config_idx_to_plot}.U_dyn,1),...
    abs(result{ismember(animals,animal)}.resdynfit_final_raw{config_idx_to_plot}.U_dyn(:,1:dynsub_dim)));
axis tight;colormap(flipud(cbrewer('div','RdYlBu',64)));
cbh = colorbar;
cbh.Location = 'eastoutside';

%% Plotting the eigenvalue magnitude/time-constants time-course (Fig 4b)

marker_size = 3;
line_width = 0.75;
err_bar_length = 0.025;
err_bar_width = 1.25;
t_idx_ci = {[0.24 0.74],[-0.3 0.05]};
face_alpha = 0;
do_color_edge = true;
markersize = 6;
ah = subplot(n_rows,n_cols,[2 3 10 11]);hold(ah);set(ah,'plotbox',[1 1 1]);
time_abs = result{ismember(animals,animal)}.resdynfit_final_raw{config_idx_to_plot}.final_model{cidx}.time;
time_rel = result{ismember(animals,animal)}.resdynfit_final_raw{config_idx_to_plot}.final_model{cidx}.time_rel;
time_labels = result{ismember(animals,animal)}.resdynfit_final_raw{config_idx_to_plot}.final_model{cidx}.time_iev;
[timeOfEvents] = findEventMarkers(time_abs,time_rel,time_labels);
col_lines = flipud(cbrewer('seq','GnBu',20));
plotpars.col_lines = col_lines(round(linspace(1, size(col_lines,1) - 10, size(result{ismember(animals,animal)}.resdynfit_final_raw{config_idx_to_plot}.final_model{cidx}.A,1))),:);
plotpars.do_plot_markers = true;
plotpars.markersize = marker_size;
plotpars.linewidth = line_width;
plotpars.ylim = [0 1.1];
plotpars.ytick = [0 : 0.2 : 1.1];
plotpars.yticklabel = [0 : 0.2 : 1.1];
plotpars.yLabel = '| eigen-value | (a.u)';
plotpars.xLabel = 'time (s)';
plotpars.ci_type = 'bar';
plotpars.do_plot_stability = true;
data_to_plot = num2cell(abs(result{ismember(animals,animal)}.resdynfit_final_raw{config_idx_to_plot}.final_model{cidx}.eigVal),2);
plotpars.line_style = repmat({'-'},[numel(data_to_plot),1]);
plotpars.marker_style = repmat('o',[length(data_to_plot) 1]);
plotpars.markerfacecol = repmat([1 1 1],[length(data_to_plot) 1]);

ci = permute(result{ismember(animals,animal)}.resdynfit_final_raw{config_idx_to_plot}.final_model{cidx}.boot_stats.eigVal,[2 3 1]);
[ci_masked] = create_bootci_masks(cell2mat(data_to_plot),ci, t_idx_ci, time_rel, time_labels);            
plotpars.bar_length = err_bar_length;
plotpars.bar_width = err_bar_width;
ci_to_plot = num2cell(permute(ci_masked,[2 3 1]),[1 2]);

[ah] = plotXvsTime(ah,data_to_plot,time_abs,plotpars,'time_iev',time_labels,'ci',ci_to_plot);
[ah] = formatXvsTimeplot(ah,timeOfEvents,plotpars.do_plot_stability);
yyaxis right
taus_to_show = [-0.045 -0.1 -0.15 -0.2 -0.35 -0.5 -1 Inf 0.5];
yticks_tau = exp((time_abs(2) - time_abs(1))./taus_to_show);
set(ah,'YColor',[0 0 0],'ylim',plotpars.ylim,'ytick',yticks_tau)
y_ticks_tau = arrayfun( @(x) sprintf('%.2f',x), taus_to_show,'uniformoutput', false);
set(ah,'yticklabel',y_ticks_tau);
set(ah,'Parent',gcf);
%% Plotting the eigenvalue phase/rotation-freq time-course (Fig 4c)

ah = subplot(n_rows,n_cols,[4 5 12 13]);hold(ah);set(ah,'plotbox',[1 1 1]);
time_abs = result{ismember(animals,animal)}.resdynfit_final_raw{config_idx_to_plot}.final_model{cidx}.time;
time_labels = result{ismember(animals,animal)}.resdynfit_final_raw{config_idx_to_plot}.final_model{cidx}.time_iev;
col_lines = flipud(cbrewer('seq','GnBu',20));
plotpars.col_lines = col_lines(round(linspace(1, size(col_lines,1) - 10, size(result{ismember(animals,animal)}.resdynfit_final_raw{config_idx_to_plot}.final_model{cidx}.A,1))),:);
plotpars.do_plot_markers = true;
plotpars.markersize = marker_size;
plotpars.linewidth = line_width;
plotpars.ylim = [-0.5 0.5];
plotpars.ytick = [-0.5 : 0.25 : 0.5];
plotpars.yticklabel = [-0.5 : 0.25 : 0.5];
plotpars.yLabel ='Angle (eigen-value) (a.u)';
plotpars.xLabel = 'time (s)';
plotpars.ci_type = 'bar';
plotpars.do_plot_stability = false;
data_to_plot = num2cell(result{ismember(animals,animal)}.resdynfit_final_raw{config_idx_to_plot}.final_model{cidx}.eigAngle,2);
plotpars.line_style = repmat({'-'},[numel(data_to_plot),1]);
plotpars.marker_style = repmat('o',[length(data_to_plot) 1]);
plotpars.markerfacecol = repmat([1 1 1],[length(data_to_plot) 1]);

ci = permute(result{ismember(animals,animal)}.resdynfit_final_raw{config_idx_to_plot}.final_model{cidx}.boot_stats.eigAngle,[2 3 1]);
[ci_masked] = create_bootci_masks(cell2mat(data_to_plot),ci, t_idx_ci, time_rel, time_labels);            
plotpars.bar_length = err_bar_length;
plotpars.bar_width = err_bar_width;
ci_to_plot = num2cell(permute(ci_masked,[2 3 1]),[1 2]);

[ah] = plotXvsTime(ah,data_to_plot,time_abs,plotpars,'time_iev',time_labels,'ci',ci_to_plot);
[ah] = formatXvsTimeplot(ah,timeOfEvents,plotpars.do_plot_stability); 
yyaxis right
set(ah,'YColor',[0 0 0],'ylim',plotpars.ylim,'ytick',plotpars.ytick)
f_to_show = [-2 -1 -0.5 -0.25 0 0.25 0.5 1 2];
yticks_f = 2*pi*(time_abs(2) - time_abs(1))*f_to_show;
set(ah,'YColor',[0 0 0],'ylim',plotpars.ylim,'ytick',yticks_f)
y_ticks_f = arrayfun( @(x) sprintf('%.2f',x), f_to_show ,'uniformoutput', false);
set(ah,'yticklabel',y_ticks_f);
set(ah,'Parent',gcf);

%% Plotting the singular value time-course (Fig 4d)

ah = subplot(n_rows,n_cols,[6 7 14 15]);hold(ah);set(ah,'plotbox',[1 1 1]);
time_abs = result{ismember(animals,animal)}.resdynfit_final_raw{config_idx_to_plot}.final_model{cidx}.time;
time_labels = result{ismember(animals,animal)}.resdynfit_final_raw{config_idx_to_plot}.final_model{cidx}.time_iev;
col_lines = flipud(cbrewer('seq','GnBu',20));
plotpars.col_lines = col_lines(round(linspace(1, size(col_lines,1) - 10, size(result{ismember(animals,animal)}.resdynfit_final_raw{config_idx_to_plot}.final_model{cidx}.A,1))),:);
plotpars.do_plot_markers = true;
plotpars.markersize = marker_size;
plotpars.linewidth = line_width;
plotpars.ylim = [0 1.1];
plotpars.ytick = [0 : 0.2 : 1.1];
plotpars.yticklabel = [0 : 0.2 : 1.1];
plotpars.yLabel = 'singular-value (a.u)';
plotpars.xLabel = 'time (s)';
plotpars.ci_type = 'bar';
plotpars.do_plot_stability = true;
data_to_plot = num2cell(abs(result{ismember(animals,animal)}.resdynfit_final_raw{config_idx_to_plot}.final_model{cidx}.sVal),2);
plotpars.line_style = repmat({'-'},[numel(data_to_plot),1]);
plotpars.marker_style = repmat('o',[length(data_to_plot) 1]);
plotpars.markerfacecol = repmat([1 1 1],[length(data_to_plot) 1]);
ci = permute(result{ismember(animals,animal)}.resdynfit_final_raw{config_idx_to_plot}.final_model{cidx}.boot_stats.sVal,[2 3 1]);
[ci_masked] = create_bootci_masks(cell2mat(data_to_plot),ci, t_idx_ci, time_rel, time_labels);            
plotpars.bar_length = err_bar_length;
plotpars.bar_width = err_bar_width;
ci_to_plot = num2cell(permute(ci_masked,[2 3 1]),[1 2]);
[ah] = plotXvsTime(ah,data_to_plot,time_abs,plotpars,'time_iev',time_labels,'ci',ci_to_plot);
[ah] = formatXvsTimeplot(ah,timeOfEvents,plotpars.do_plot_stability);

%% Plotting the ev/sv summary for all animals (Fig 4e)

subplot_idxs1 = {[25 26 33 34];[41 42 49 50]};
subplot_idxs2 = {[27 28 35 36];[43 44 51 52]};
subplot_idxs3 = {[29 30 37 38];[45 46 53 54]};
subplot_idxs4 = {[31 32 39 40];[47 48 55 56]};
marker_size = 3;
line_width = 0.75;

for iani = 1: length(animals)
    
    % Plotting maximum ev magnitude/time-constants
    clear plotpars
    n_choices = 2;
    n_configs = length(result{iani}.resdynfit_final_raw);
    ah = subplot(n_rows,n_cols,subplot_idxs1{iani});hold(ah);set(ah,'plotbox',[1 1 1]);
    time_abs = result{iani}.resdynfit_final_raw{config_idx_to_plot}.final_model{cidx}.time;
    time_rel = result{iani}.resdynfit_final_raw{config_idx_to_plot}.final_model{cidx}.time_rel;
    time_labels = result{iani}.resdynfit_final_raw{config_idx_to_plot}.final_model{cidx}.time_iev;
    [timeOfEvents] = findEventMarkers(time_abs,time_rel,time_labels);
    col_lines_pref = [0 0 1];
    col_lines_anti = [1 0 0];
    config_markers = {'^';'d';'o';'s'};
    plotpars.do_plot_markers = true;
    plotpars.markersize = marker_size;
    plotpars.linewidth = line_width;
    plotpars.ylim = [0 1.1];
    plotpars.ytick = [0 : 0.2 : 1.1];
    plotpars.yticklabel = [0 : 0.2 : 1.1];
    plotpars.yLabel = '| max eigen-value | (a.u)';
    plotpars.xLabel = 'time (s)';
    plotpars.do_plot_stability = true;
    data_to_plot = cell(n_choices*n_configs,1);
    c_count = 1;
    for icond = 1:n_choices
        for iconfig = 1:n_configs

            if(iconfig == 1)
               if(icond == 1)
                  plotpars.col_lines(c_count,:) = col_lines_anti; 

               else
                   plotpars.col_lines(c_count,:) = col_lines_pref; 

               end

            else
                if(icond == 1)
                    plotpars.col_lines(c_count,:) = col_lines_pref;

                else
                    plotpars.col_lines(c_count,:) = col_lines_anti;

                end

            end
            data_to_plot{c_count} = max(result{iani}.resdynfit_final_raw{iconfig}.final_model{icond}.eigVal,[],1);
            plotpars.marker_style(c_count) = config_markers{iconfig};

            c_count = c_count + 1;
        end

    end
    data_to_plot = cat(1,data_to_plot,{mean(cell2mat(data_to_plot),1)});
    plotpars.col_lines = cat(1,plotpars.col_lines,[0 0 0]);
    plotpars.marker_style = cat(1,plotpars.marker_style',' ');
    plotpars.line_style = repmat({'-'},[numel(data_to_plot),1]);
    plotpars.markerfacecol = repmat([1 1 1],[length(data_to_plot) 1]);
    [ah] = plotXvsTime(ah,data_to_plot,time_abs,plotpars,'time_iev',time_labels);
    [ah] = formatXvsTimeplot(ah,timeOfEvents,plotpars.do_plot_stability);
    yyaxis right
    taus_to_show = [-0.045 -0.1 -0.15 -0.2 -0.35 -0.5 -1 Inf 0.5];
    yticks_tau = exp((time_abs(2) - time_abs(1))./taus_to_show);
    set(ah,'YColor',[0 0 0],'ylim',plotpars.ylim,'ytick',yticks_tau)
    y_ticks_tau = arrayfun( @(x) sprintf('%.2f',x), taus_to_show,'uniformoutput', false);
    set(ah,'yticklabel',y_ticks_tau);
    set(ah,'Parent',gcf);
    
    % Plotting maximum ev phase/rot-freq
    clear plotpars
    n_choices = 2;
    ah = subplot(n_rows,n_cols,subplot_idxs2{iani});hold(ah);set(ah,'plotbox',[1 1 1]);
    time_abs = result{iani}.resdynfit_final_raw{config_idx_to_plot}.final_model{cidx}.time;
    time_labels = result{iani}.resdynfit_final_raw{config_idx_to_plot}.final_model{cidx}.time_iev;
    col_lines_pref = [0 0 1];
    col_lines_anti = [1 0 0];
    config_markers = {'^';'d';'o';'s'};
    plotpars.do_plot_markers = true;
    plotpars.markersize = marker_size;
    plotpars.linewidth = line_width;
    plotpars.ylim = [0 0.8];
    plotpars.ytick = [0: 0.2 : 0.8];
    plotpars.yticklabel = [0: 0.2 : 0.8];
    plotpars.yLabel = 'max angle (eigen-value) (a.u)';
    plotpars.xLabel = 'time (s)';
    plotpars.do_plot_stability = false;
    data_to_plot = cell(n_choices*n_configs,1);
    c_count = 1;
    for icond = 1:n_choices
        for iconfig = 1:n_configs

            if(iconfig == 1)
               if(icond == 1)
                  plotpars.col_lines(c_count,:) = col_lines_anti; 

               else
                   plotpars.col_lines(c_count,:) = col_lines_pref; 

               end

            else
                if(icond == 1)
                    plotpars.col_lines(c_count,:) = col_lines_pref;

                else
                    plotpars.col_lines(c_count,:) = col_lines_anti;

                end

            end
            data_to_plot{c_count} = max(abs(result{iani}.resdynfit_final_raw{iconfig}.final_model{icond}.eigAngle),[],1);
            plotpars.marker_style(c_count) = config_markers{iconfig};

            c_count = c_count + 1;
        end

    end
    data_to_plot = cat(1,data_to_plot,{mean(cell2mat(data_to_plot),1)});
    plotpars.col_lines = cat(1,plotpars.col_lines,[0 0 0]);
    plotpars.marker_style = cat(1,plotpars.marker_style',' ');
    plotpars.line_style = repmat({'-'},[numel(data_to_plot),1]);
    plotpars.markerfacecol = repmat([1 1 1],[length(data_to_plot) 1]);
    [ah] = plotXvsTime(ah,data_to_plot,time_abs,plotpars,'time_iev',time_labels);
    [ah] = formatXvsTimeplot(ah,timeOfEvents,plotpars.do_plot_stability);
    yyaxis right
    set(ah,'YColor',[0 0 0],'ylim',plotpars.ylim,'ytick',plotpars.ytick)
    f_to_show = [0 0.25 0.5 1 2 2.5];
    yticks_f = 2*pi*(time_abs(2) - time_abs(1))*f_to_show;
    set(ah,'YColor',[0 0 0],'ylim',plotpars.ylim,'ytick',yticks_f)
    y_ticks_f = arrayfun( @(x) sprintf('%.2f',x), f_to_show ,'uniformoutput', false);
    set(ah,'yticklabel',y_ticks_f);
    set(ah,'Parent',gcf);
    
    % Plotting the maximum singular value
    clear plotpars
    n_choices = 2;
    ah = subplot(n_rows,n_cols,subplot_idxs3{iani});hold(ah);set(ah,'plotbox',[1 1 1]);
    time_abs = result{iani}.resdynfit_final_raw{config_idx_to_plot}.final_model{cidx}.time;
    time_labels = result{iani}.resdynfit_final_raw{config_idx_to_plot}.final_model{cidx}.time_iev;
    col_lines_pref = [0 0 1];
    col_lines_anti = [1 0 0];
    config_markers = {'^';'d';'o';'s'};
    plotpars.do_plot_markers = true;
    plotpars.markersize = marker_size;
    plotpars.linewidth = line_width;
    plotpars.ylim = [0 1.1];
    plotpars.ytick = [0 : 0.2 : 1.1];
    plotpars.yticklabel = [0 : 0.2 : 1.1];
    plotpars.yLabel = 'max singular value (a.u)';
    plotpars.xLabel = 'time (s)';
    plotpars.do_plot_stability = true;
    data_to_plot = cell(n_choices*n_configs,1);
    c_count = 1;
    for icond = 1:n_choices
        for iconfig = 1:n_configs

            if(iconfig == 1)
               if(icond == 1)
                  plotpars.col_lines(c_count,:) = col_lines_anti; 

               else
                   plotpars.col_lines(c_count,:) = col_lines_pref; 

               end

            else
                if(icond == 1)
                    plotpars.col_lines(c_count,:) = col_lines_pref;

                else
                    plotpars.col_lines(c_count,:) = col_lines_anti;

                end

            end
            data_to_plot{c_count} = max(result{iani}.resdynfit_final_raw{iconfig}.final_model{icond}.sVal,[],1);
            plotpars.marker_style(c_count) = config_markers{iconfig};

            c_count = c_count + 1;
        end

    end
    data_to_plot = cat(1,data_to_plot,{mean(cell2mat(data_to_plot),1)});
    plotpars.col_lines = cat(1,plotpars.col_lines,[0 0 0]);
    plotpars.marker_style = cat(1,plotpars.marker_style',' ');
    plotpars.line_style = repmat({'-'},[numel(data_to_plot),1]);
    plotpars.markerfacecol = repmat([1 1 1],[length(data_to_plot) 1]);
    [ah] = plotXvsTime(ah,data_to_plot,time_abs,plotpars,'time_iev',time_labels);
    [ah] = formatXvsTimeplot(ah,timeOfEvents,plotpars.do_plot_stability);
    
    % Plotting the non-normality measure
    clear plotpars
    n_choices = 2;
    ah = subplot(n_rows,n_cols,subplot_idxs4{iani});hold(ah);set(ah,'plotbox',[1 1 1]);
    time_abs = result{iani}.resdynfit_final_raw{config_idx_to_plot}.final_model{cidx}.time;
    time_labels = result{iani}.resdynfit_final_raw{config_idx_to_plot}.final_model{cidx}.time_iev;
    col_lines_pref = [0 0 1];
    col_lines_anti = [1 0 0];
    config_markers = {'^';'d';'o';'s'};
    plotpars.do_plot_markers = true;
    plotpars.markersize = marker_size;
    plotpars.linewidth = line_width;
    plotpars.ylim = [0 0.5];
    plotpars.ytick = [0 : 0.2 : 0.5];
    plotpars.yticklabel = [0 : 0.2 : 0.5];
    plotpars.yLabel = 'non-normality (a.u)';
    plotpars.xLabel = 'time (s)';
    plotpars.do_plot_stability = true;
    data_to_plot = cell(n_choices*n_configs,1);
    c_count = 1;
    for icond = 1:n_choices
        for iconfig = 1:n_configs

            if(iconfig == 1)
               if(icond == 1)
                  plotpars.col_lines(c_count,:) = col_lines_anti; 

               else
                   plotpars.col_lines(c_count,:) = col_lines_pref; 

               end

            else
                if(icond == 1)
                    plotpars.col_lines(c_count,:) = col_lines_pref;

                else
                    plotpars.col_lines(c_count,:) = col_lines_anti;

                end

            end
            data_to_plot{c_count} = result{iani}.resdynfit_final_raw{iconfig}.final_model{icond}.depN;
            plotpars.marker_style(c_count) = config_markers{iconfig};

            c_count = c_count + 1;
        end

    end
    data_to_plot = cat(1,data_to_plot,{mean(cell2mat(data_to_plot),1)});
    plotpars.col_lines = cat(1,plotpars.col_lines,[0 0 0]);
    plotpars.marker_style = cat(1,plotpars.marker_style',' ');
    plotpars.line_style = repmat({'-'},[numel(data_to_plot),1]);
    plotpars.markerfacecol = repmat([1 1 1],[length(data_to_plot) 1]);
    [ah] = plotXvsTime(ah,data_to_plot,time_abs,plotpars,'time_iev',time_labels);
    [ah] = formatXvsTimeplot(ah,timeOfEvents,plotpars.do_plot_stability);
    
end
