%{
This script fits the linear models of the following form:

angle_with_x = b0 + b1*|ev| + b2*(rot_freq)
where, x can either be choice, time, JPC12 or jPC34.

The fits are done separately for each epoch (stim or sacc), by pooling
results acrosss choices, dimensions and task configurations.

It plots the estimated coefficients for the ev (b1) and rot_freq (b2), and
also plots the raw data (angle_with_x vs ev/rot_freq) and the corresponding
fits.

See Figure 5
%}

%% Clear workspace

clearvars -except DIRS
close all
clc
%% Load data

animal = 'Tex';
bin_size = 45;
cd('/Users/Aniruddh/Work_PhD/Residual Dynamics/residual_dynamics');
data_path = './data/processedneuraldata';
data_file_name = sprintf('%s%s',animal,'_ndim=8_lag=3_alpha=200_50_allconfigsresdynresults.mat');
load(fullfile(data_path, data_file_name))

% Properties
dt =  bin_size/1000;
save_fig_path = ''; % enter a valid path;
do_save_fig = false;
do_mean_subtract = true;
do_standardize = false;

%% AUGMENTING DATA-TABLE WITH AUXILLARY VARIABLES

data_tbl = result.resdynfit;
alignments = result.cond_avg_pars.jPC.alignments;
data_tbl.abs_rot_freq = abs(data_tbl.rot_freq); % store the absolute value of rotation frequency
unique_time_labels = unique(data_tbl.time_label);
n_alignments = length(unique_time_labels);

%% Defining and fitting all the linear models

lm_model_specs = {'1', '1 + ev', '1 + abs_rot_freq',...
    '1 + ev + abs_rot_freq', '1 + ev + abs_rot_freq + ev:abs_rot_freq'};
lm_model_names = {'n','n + ev','n + rf','n + ev + rf','n + ev + rf + ev*rf'};
response_variables = {'angle_with_choice','angle_with_time', 'angle_with_jPC12', 'angle_with_jPC34'};
planes_names = cellfun(@(x) strrep(x,'angle_with_',''),response_variables,'uni',false);
n_models = length(lm_model_names);
lm_mdls = cell(n_alignments, length(response_variables));
means = cell(n_alignments,1); % when mean-centering, you need to store the means that are subtracted

for iev = 1: n_alignments
    data_tbl_ev = data_tbl(data_tbl.time_label == iev,:);
    
    if(do_mean_subtract)
        means{iev}.ev = mean(data_tbl_ev.ev);
        means{iev}.abs_rot_freq = mean(data_tbl_ev.abs_rot_freq);
        data_tbl_ev.ev = data_tbl_ev.ev - means{iev}.ev;
        data_tbl_ev.abs_rot_freq = data_tbl_ev.abs_rot_freq - means{iev}.abs_rot_freq;
         
    end
    
    if(do_standardize)
        data_tbl_ev.ev = zscore(data_tbl_ev.ev);
        data_tbl_ev.abs_rot_freq = zscore(data_tbl_ev.abs_rot_freq);
    end
    
    for i_rv = 1: length(response_variables)
        lm_mdls{iev, i_rv} = cell(n_models,1);
        for imdl = 1 : n_models
       
            model_formula = sprintf('%s%s%s',response_variables{i_rv} , '~', lm_model_specs{imdl});     
            lm_mdls{iev, i_rv}{imdl} = fitlm(data_tbl_ev, model_formula);
                
        end
        
    end
    
end


%% Plotting the coefficients of the linear model fits and the scatter plots

clear plotpars
model_type_to_plot = '1 + ev + abs_rot_freq';
idx_model_simple = find(ismember(lm_model_specs,model_type_to_plot));

mdls_to_plot = cellfun(@(x) x{idx_model_simple},lm_mdls,'uni',false);
plotpars.coef_plot.coef_cols = [0 0 0; 0 0 1];
plotpars.coef_plot.coefs_to_plot = {'ev','abs_rot_freq'};
plotpars.plane_names = planes_names;
plotpars.alignment_labels = alignments;
plotpars.fig_identifier = sprintf('%s%s%d%s',animal,'_binsize=',bin_size,'ms');

figure; ah = gca; hold(ah); set(ah,'plotbox',[1 1 1]); 
plot_coefficients(ah, mdls_to_plot, plotpars)
if(do_save_fig)
    figname = sprintf('%s%s','LinModelCoefs_',data_file_name);
    export_fig(fullfile(save_fig_path,figname),'-painters','-transparent','-pdf');
end

xvars_to_plot = plotpars.coef_plot.coefs_to_plot;
plotpars.scatter.marker_size = 4;
plotpars.scatter.marker_col = [0.5 0.5 0.5];
plotpars.scatter.mdl_fit_style = '--';

for iev = 1: n_alignments
    h = figure;set(h,'Position',[102 109 1026 585]);
    data_tbl_ev = data_tbl(data_tbl.time_label == iev,:);
    
    if(do_mean_subtract)
        
        data_tbl_ev.ev = data_tbl_ev.ev - means{iev}.ev;
        data_tbl_ev.abs_rot_freq = data_tbl_ev.abs_rot_freq - means{iev}.abs_rot_freq;
         
    end
    
    for ivar = 1: length(xvars_to_plot)
      for i_rv = 1: length(response_variables)
          
          ah = subplot(length(xvars_to_plot),length(response_variables),(ivar-1)*length(response_variables)+i_rv);
          hold(ah);set(ah,'plotbox',[1 1 1]);
          [ah] = plot_fit_and_scatter(ah, data_tbl_ev.(xvars_to_plot{ivar}), data_tbl_ev.(response_variables{i_rv}),...
              xvars_to_plot{ivar}, mdls_to_plot{iev,i_rv}, plotpars, means{iev}.(xvars_to_plot{ivar}));

          ylabel(ah, response_variables{i_rv},'Interpreter','none')                 
          set(ah,'Parent',h);
          
      end
    end
    
    suptitle(sprintf('%s%s%s%s',animal, ' ,', plotpars.alignment_labels{iev}, ' epoch'))
    if(do_save_fig)
        figname = sprintf('%s%s%s%s','ScatterandfitsLinModel_', plotpars.alignment_labels{iev},'_',data_file_name);
        export_fig(fullfile(plotpars.save_fig_path,figname),'-painters','-transparent','-pdf');
    end
end
%% Computing and plotting the goodness of fits for the the different linear models

BIC = cell(n_alignments, length(response_variables));
Rsquared = cell(n_alignments, length(response_variables));

for iev = 1: n_alignments
    for i_rv = 1: length(response_variables)
        
        BIC{iev,i_rv} = cell2mat(cellfun(@(x) x.ModelCriterion.BIC,lm_mdls{iev,i_rv},'uni',false));
        Rsquared{iev,i_rv} = cell2mat(cellfun(@(x) x.Rsquared.Ordinary,lm_mdls{iev,i_rv},'uni',false));
        
    end
end

% Plotting BIC
for iev = 1: n_alignments
    figure; 
    for i_rv = 1: length(response_variables)
        
        
        ah = subplot(1,length(response_variables),i_rv);hold(ah);set(ah,'plotbox',[1 4 1]);
        set(ah,'ylim',[0.5 n_models + 0.5],'ytick',1:n_models);
        ylims = get(ah,'ylim');
    

        plot(ah,[0 0],ylims,'--','Color',[0.6 0.6 0.6]);
        
        delta_bic = BIC{iev, i_rv} -  BIC{iev, i_rv}(1);
        [min_del_bic,idx] = min(delta_bic);
        
        
        plot(ah,[min_del_bic-10 min_del_bic-10],ylims,'--','Color',[1 0 0]);
        plot(ah,[min_del_bic+10 min_del_bic+10],ylims,'--','Color',[1 0 0]);
        
        if(i_rv == 1)
            set(ah,'yticklabel',lm_model_names);
            ylabel(ah, '(model types)')
        else
            set(ah,'yticklabel',{});
        end

        plot(ah,BIC{iev, i_rv} -  BIC{iev, i_rv}(1),1:n_models, 'o','markerfacecolor',[ 0 0 0],'markeredgecolor',[0 0 0],'markersize',5);
   
        title(ah,planes_names{i_rv});
        xlabel(ah, '\Delta BIC');
    end
    
    if(do_save_fig)
        fig_name = sprintf('%s%s%s%s','deltaBIC_',alignments{iev},'_',strrep(file_name,'.mat',''));
        export_fig(fullfile(save_fig_path,fig_name),'-painters','-transparent','-pdf');

    end
    
end


% Plotting Rsquared

for iev = 1: n_alignments   

    figure; 
    for i_rv = 1: length(response_variables)
        
        
        ah = subplot(1,length(response_variables),i_rv);hold(ah);set(ah,'plotbox',[1 4 1]);
        set(ah,'ylim',[0.5 n_models + 0.5],'ytick',1:n_models);
        ylims = get(ah,'ylim');        
        if(i_rv == 1)
            set(ah,'yticklabel',lm_model_names);
            ylabel(ah, '(model types)')
        else
            set(ah,'yticklabel',{});
        end

        plot(ah,Rsquared{iev, i_rv},1:n_models, 'o','markerfacecolor',[ 0 0 0],'markeredgecolor',[0 0 0],'markersize',5);
        xlabel(ah, 'R^{2}');

        title(ah,planes_names{i_rv});

    end
    
    if(do_save_fig)
        fig_name = sprintf('%s%s%s%s','Rsquared_',alignments{iev},'_',strrep(file_name,'.mat',''));
        export_fig(fullfile(save_fig_path,fig_name),'-painters','-transparent','-pdf');
    end
end