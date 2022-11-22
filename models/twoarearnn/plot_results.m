clearvars -except DIRS
clc
close all

%% Loading stored results Select the networks and conditions that you want to plot

% Specifying which models and choice goes into the final plot
% This specifies all the relevant model configurations for both types of
% networks (i.e with or without feedback)
model_pars = table([0.3600;0.2000],...
                   [0.0800;0.2400],...
                   [0;0.2400],...
             'VariableNames',{'self' 'fwd' 'fbk'});

file_path = './data/twoarearnn/';
file_name = 'nofeedbackmodels.mat'; % choose either 'withfeedback' or 'nofeedback' here
do_save_vars = false;

[result,self_rec,fwd,fbk, model_out] = extract_model_analyses(file_path,file_name,model_pars);
cidx = 1;  % Choice condition to plot (plots choice-1 when cidx = 1)
degenerate_idxs = (fwd == 0 & fbk == 0); % indices which correspond to degenerate models

%% Defining the task-relevant dimensions of the model in observation space

task_rel_axes_obs.ppc.choice  = result.obs_pars.C_obs*normc([1;-1;0;0]);
task_rel_axes_obs.pfc.choice  = result.obs_pars.C_obs*normc([0;0;1;-1]);
task_rel_axes_obs.ppc.time = result.obs_pars.C_obs*normc([1;1;0;0]);
task_rel_axes_obs.pfc.time = result.obs_pars.C_obs*normc([0;0;1;1]);

%% Plotting results

% colors for eigenvalue plots
col_lines = [8 64 129; 9 109 175; 48 145 193; 60 180 223]./255;

% colors for task overlap plots
cols_ppc = [0.0430    0.3555    0.1016;...
    0.4609    0.1797    0.6836];
cols_hsv = rgb2hsv(cols_ppc);
cols_hsv(1, 2) = cols_hsv(1, 2) .* 0.1;
cols_hsv(1, 1) = cols_hsv(1, 1) .* 1.15;
cols_hsv(2, 2) = cols_hsv(2, 2) .* 0.4;
cols_pfc = hsv2rgb(cols_hsv);
cols = [cols_ppc;cols_pfc];

ms_size= 10; %marker size

% properites of error bars
error_bar_length = 0.01;
error_bar_width = 0.5;
error_line_width = 0.5;

% Additional properties for computing overlap with task axes
vars_to_store = {'angle'};
area_to_use_for_tmax = 'pfc';
dim_masks = {[true(result.sim_pars.obs_dim,1); false(result.sim_pars.obs_dim,1)];...
    [false(result.sim_pars.obs_dim,1); true(result.sim_pars.obs_dim,1)];
    [true(result.sim_pars.obs_dim,1) ; true(result.sim_pars.obs_dim,1)]};

% Making all the relevant plots for the example models specified in
% model_pars
for imdl = 1: length(model_out)
    
        h1 = figure;
        h2= figure;
        h3 = figure;
        % Plot the condition-averaged trajectory
        
        % ppc
        ah = subplot(1,2,1);hold(ah);set(ah,'plotbox',[1,1 1],'Parent',h1);
        plot(ah,model_out{imdl}.data.task_modes.ppc.choice.proj(1,:),model_out{imdl}.data.task_modes.ppc.time.proj(1,:),'-','Color',[0 0 1]);
        plot(ah,model_out{imdl}.data.task_modes.ppc.choice.proj(2,:),model_out{imdl}.data.task_modes.ppc.time.proj(2,:),'-','Color',[1 0 0]);
        plot(ah,model_out{imdl}.data.task_modes.ppc.choice.proj(1,:),model_out{imdl}.data.task_modes.ppc.time.proj(1,:),'o','markeredgecolor',[0 0 1],'markerfacecolor',[1 1 1]);
        plot(ah,model_out{imdl}.data.task_modes.ppc.choice.proj(2,:),model_out{imdl}.data.task_modes.ppc.time.proj(2,:),'o','markeredgecolor',[1 0 0],'markerfacecolor',[1 1 1]);
        xlabel(ah,'choice-mode (a.u)')
        ylabel(ah,'time-mode (a.u)')
        
        % pfc
        ah = subplot(1,2,2);hold(ah);set(ah,'plotbox',[1,1 1],'Parent',h1);
        plot(ah,model_out{imdl}.data.task_modes.pfc.choice.proj(1,:),model_out{imdl}.data.task_modes.pfc.time.proj(1,:),'-','Color',[0 0 1]);
        plot(ah,model_out{imdl}.data.task_modes.pfc.choice.proj(2,:),model_out{imdl}.data.task_modes.pfc.time.proj(2,:),'-','Color',[1 0 0]);
        plot(ah,model_out{imdl}.data.task_modes.pfc.choice.proj(1,:),model_out{imdl}.data.task_modes.pfc.time.proj(1,:),'o','markeredgecolor',[0 0 1],'markerfacecolor',[1 1 1]);
        plot(ah,model_out{imdl}.data.task_modes.pfc.choice.proj(2,:),model_out{imdl}.data.task_modes.pfc.time.proj(2,:),'o','markeredgecolor',[1 0 0],'markerfacecolor',[1 1 1]);
        xlabel(ah,'choice-mode (a.u)')
        ylabel(ah,'time-mode (a.u)')
        
        % Plot the local and global residual dynamics (eigenvalue only) of the model
        for iarea = 1: length(model_out{imdl}.result_cv)

            ah = subplot(1, 3, iarea);
            hold(ah); set(ah,'plotbox',[1 1 1])
            
            time_abs = model_out{imdl}.result_cv{iarea}.resdyn.final_model{cidx}.time_rel;
            time_labels = model_out{imdl}.result_cv{iarea}.resdyn.final_model{cidx}.time_iev;
            clear plotpars
            plotpars.col_lines = col_lines;
            plotpars.do_plot_markers = true;
            plotpars.markersize = 3;
            plotpars.linewidth = 0.75;
            plotpars.ylim = [0 1.4];
            plotpars.ytick = [0 : 0.25 : 1.25];
            plotpars.yticklabel = [0 : 0.25 : 1.25];
            plotpars.yLabel = '| eigen-value | (a.u)';
            plotpars.xLabel = 'time (s)';
            plotpars.ci_type = 'bar';
            plotpars.do_plot_stability = true;
            data_to_plot = num2cell(abs(model_out{imdl}.result_cv{iarea}.resdyn.final_model{cidx}.eigVal),2);
            plotpars.line_style = repmat({'-'},[numel(data_to_plot),1]);
            plotpars.marker_style = repmat('o',[length(data_to_plot) 1]);
            plotpars.markerfacecol = repmat([1 1 1],[length(data_to_plot) 1]);
            [ah] = plotXvsTime(ah,data_to_plot,time_abs,plotpars,'time_iev',time_labels);
            [ah] = formatXvsTimeplot(ah,[],plotpars.do_plot_stability);
            title(ah,model_out{imdl}.areas{iarea})
            set(ah,'Parent',h2);
        end
        
        % compute overlap of eigenvectors with task relevant dimensions
        % defined separately for each area
        
        task_overlap = struct();
        for iarea = 1: length(model_out{imdl}.result_cv)
            num_dims =  size(model_out{imdl}.result_cv{iarea}.resdyn.final_model{cidx}.A,1);
            U_dyn = model_out{imdl}.result_cv{iarea}.resdyn.U_dyn(:,1:num_dims);
            
            for tt = 1:length(model_out{imdl}.result_cv{iarea}.resdyn.final_model{cidx}.time_rel)
                eval = model_out{imdl}.result_cv{iarea}.resdyn.final_model{cidx}.eigVal(:,tt);
                evec = model_out{imdl}.result_cv{iarea}.resdyn.final_model{cidx}.eigVec(:,:,tt);
                
                [task_overlap.(model_out{imdl}.areas{iarea}){tt}] = ...
                    compute_eigenvector_taskrelaxes_overlap(evec, eval, U_dyn,...
                    task_rel_axes_obs, dim_masks{iarea}, vars_to_store);
                    
            end
                        
        end
        
        % Evaluting the time at which the ev magnitude of the local dynamics 
        % in the area specified by the variable 'area_to_use_for_tmax'
        % reaches its maximum
        idx_area_tmax = ismember(model_out{imdl}.areas,area_to_use_for_tmax);
        tmax = extract_timeidx_of_maxev(model_out{imdl}.result_cv{idx_area_tmax}.resdyn);
        
        
        for iarea = 1: length(model_out{imdl}.result_cv)    
            
            % Selecting all angles to plot for the global resdyn 
            
            if(strcmp(model_out{imdl}.areas{iarea},'both'))

                angles_to_plot = {task_overlap.(model_out{imdl}.areas{iarea}){tmax(cidx)}.ppc.choice.angle,...
                                  task_overlap.(model_out{imdl}.areas{iarea}){tmax(cidx)}.ppc.time.angle,...
                                  task_overlap.(model_out{imdl}.areas{iarea}){tmax(cidx)}.pfc.choice.angle,...
                                  task_overlap.(model_out{imdl}.areas{iarea}){tmax(cidx)}.pfc.time.angle};
            
            % Selecting all angles to plot for the local resdyn - plotting
            % angle of ppc choice/time with eigenvectors of pfc local
            % dynamics does not make sense
            elseif(strcmp(result.areas{iarea},'ppc'))
                
                angles_to_plot = {task_overlap.(model_out{imdl}.areas{iarea}){tmax(cidx)}.ppc.choice.angle,...
                  task_overlap.(model_out{imdl}.areas{iarea}){tmax(cidx)}.ppc.time.angle};
              
            elseif(strcmp(result.areas{iarea},'pfc'))
                
                angles_to_plot = {task_overlap.(model_out{imdl}.areas{iarea}){tmax(cidx)}.pfc.choice.angle,...
                  task_overlap.(model_out{imdl}.areas{iarea}){tmax(cidx)}.pfc.time.angle};
                
            end
            
            % Plotting the overlap of the task modes with the ev at tmax
            ev_at_tmax = model_out{imdl}.result_cv{iarea}.resdyn.final_model{cidx}.eigVal(:,tmax(cidx));
            [~,sorted_ev_inds] = sort(ev_at_tmax,'descend');
            
            for iev = 1:length(sorted_ev_inds)
                if(strcmp(model_out{imdl}.areas{iarea},'both'))
                    ah = subplot(2,4, iev);hold(ah);set(ah,'plotbox',[1,1 1]);
                    set(ah,'xlim',[0.5 4.5],'xtick',[],'xticklabel',[]);
                elseif(strcmp(model_out{imdl}.areas{iarea},'ppc'))
                    ah = subplot(2,4, 4 + iev);hold(ah);set(ah,'plotbox',[1,1 1]);
                    set(ah,'xlim',[0.5 2.5],'xtick',[],'xticklabel',[]);
                elseif(strcmp(model_out{imdl}.areas{iarea},'pfc'))
                    ah = subplot(2,4, 6 + iev);hold(ah);set(ah,'plotbox',[1,1 1]);
                    set(ah,'xlim',[0.5 2.5],'xtick',[],'xticklabel',[]);
                end
                
                set(ah,'ylim',[0 90],'ytick',[0 45 90],'yticklabel',[90 45 0]);
                ang_plot = cell2mat(cellfun(@(x) x(sorted_ev_inds(iev)), angles_to_plot,'uni',false));
                
                for iang = 1: length(ang_plot)
                    if(strcmp(model_out{imdl}.areas{iarea},'ppc') | strcmp(model_out{imdl}.areas{iarea},'both'))
                        stem(iang, 90 - ang_plot(iang), 'LineStyle','-', 'Color',[0 0 0],...
                            'MarkerFaceColor', cols(iang,:), 'MarkerEdgeColor',[0 0 0],'markersize',ms_size);
                    else
                        stem(iang, 90 - ang_plot(iang), 'LineStyle','-', 'Color',[0 0 0],...
                            'MarkerFaceColor', cols(2+iang,:), 'MarkerEdgeColor',[0 0 0],'markersize',ms_size);
                    end
                    
                end
                ylabel(ah,'angle(deg)')
                if(strcmp(model_out{imdl}.areas{iarea},'both'))
                    title(ah,sprintf('%s%d%s','EV',iev,', global'),'FontSize',8)
                else
                    title(ah,sprintf('%s%d%s%s','EV',iev,', local-',model_out{imdl}.areas{iarea}),'FontSize',8)
                end
                set(ah,'Parent',h3)
            end
                
        end

end

% Plotting the maximum eigenvalue along the choice mode for local dynamics
% in each area.
h4 = figure;
unique_self_vals = unique(self_rec);
count = 1;
col_mat = flipud(cbrewer('seq','Greys',8));
for iself = 1:length(unique_self_vals) - 1
    
    idx_model_self = (self_rec == unique_self_vals(iself));
    idx_model_self = idx_model_self & ~degenerate_idxs; % Filtering out the degenerate models 
    [r,c] = find(idx_model_self);
    area_names = result.areas(~ismember(result.areas,'both')); % Only plotting the local dyn.
    
    for iarea = 1:length(area_names)
        
        ah = subplot(1,2,iarea);hold all;set(ah,'plotbox',[1 1 1]);
        
        set(ah,'ylim',[0.4 2.1],'ytick',[0.5:0.5:2.0],'yticklabel',[0.5:0.5:2.0],...
            'xlim',[0 0.25],'xtick',[0:0.05:0.25],'xticklabel',[0:0.05:0.25] );
 
        max_ev_plot = NaN(length(r),1);
        
        for idx = 1: length(r)
            
            tmax = extract_timeidx_of_maxev(result.analysis{r(idx),c(idx),iarea}.resdyn_along_choice);
            
            max_ev_plot(idx) =  result.analysis{r(idx),c(idx),iarea}.resdyn_along_choice.final_model{cidx}.eigVal(tmax(cidx));
            
            n_resamples = length(result.analysis{r(idx),c(idx),iarea}.resdyn_along_choice.final_model{cidx}.boot_stats);
            max_ev_plot_boot = NaN(1,n_resamples);
            for i_resample = 1: n_resamples
                max_ev_plot_boot(i_resample) = squeeze(result.analysis{r(idx),c(idx),iarea}.resdyn_along_choice.final_model{cidx}...
                    .boot_stats{i_resample}.eigVal(tmax(cidx)));
            end
            
            plot(ah,fwd(r(idx),c(idx)), max_ev_plot(idx),'o','markerfacecolor',col_mat(count,:),'markeredgecolor',col_mat(count,:));
            
            yy = median(max_ev_plot_boot);
            ypos = prctile(max_ev_plot_boot,97.5); % Computing 95th percentiles
            yneg = prctile(max_ev_plot_boot,2.5);
            
            h_l = ploterr(fwd(r(idx),c(idx)),yy,[],{yneg ypos},'-','abshhy', error_bar_length);
            set(h_l(1),'Color',col_mat(count,:),'linewidth',error_line_width);
            set(h_l(2),'Color',col_mat(count,:),'linewidth',error_bar_width);
            
        end
        
        plot(ah,fwd(idx_model_self),max_ev_plot,'-','Color',col_mat(count,:));
        
        xlims = get(ah,'xlim');
        plot([xlims(1) xlims(2)],[1 1],'--','Color',[1 0 0]);
        title(strcat(area_names{iarea}))
        xlabel('between area connection strength')
        ylabel('eigenvalue along choice (tmax)')
    end
    count = count+1;  
end