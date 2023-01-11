%{ This script generates the different categories of phase portraits 
% (recurrent, effective and residual) for the models of integraion and
% movement (see Extended Data Fig 1 and Fig 2a-b)
%}

clearvars -except DIRS
close all 
clc
do_save_fig = false;
data_load_path = fullfile(DIRS.analysis,'/simulations/toymodels/');
save_fig_path = ''; % enter a valid path indicating where to save figures

file_name = 'integration_models_newAnalysisPipeline_09-Dec-2018.mat';
integration_models = load(fullfile(data_load_path,file_name));
file_name = 'movement_models_newAnalysisPipeline_22-Nov-2018.mat';
movement_models = load(fullfile(data_load_path,file_name));

%%
clear model_analyses
model_analyses{1} = integration_models.model_analyses;
model_analyses{2} = movement_models.model_analyses;
clear integration_models movement_models
rseed = 543612;
%%
close all;
col_markers = [0 0 1;1 0 0];
pars.box.halfWidthX = 0.25;
pars.box.halfWidthY = 0.25;
pars.time_steps = 2;
pars.box.step = 0.01;
pars.rseed = rseed;
% scale_facs = {[6 8 2.5]; [30 15 2.5]}; % original scale 
scale_facs = {[12 40 5]; [75 15 5]};  % large scale
plot_locs = {[false;false;false];[false;false;false]};
pars.arrowShape = [1 0.6 0.004 0.06 0.01 0.025 0.000 0.005]; % [k_w k_j l1 l2 a1 a2 s1 s2]
model_types = {'integration';'movement'};
t_idx = 21;
c_idx = 1;
pars.col_marker = col_markers(c_idx,:);
model_names = {{'saddle';'lineAtt';'movingFP'};{'rotational';'funnel';'movingFP'}};

for i_model = 1 : length(model_types)
    for idx_mdl = 1 : 3
        
        pars.scale_fac = scale_facs{i_model}(idx_mdl);
        pars.plot_locs = plot_locs{i_model}(idx_mdl);
        plotAllFlows(model_analyses{i_model}, model_types{i_model}, idx_mdl, c_idx, t_idx , pars)
        set(gcf,'Position',[364 851 1308 143]);
        figName = sprintf('%s%s%s%s%s','flowDiagrams_',model_types{i_model},'_',model_names{i_model}{idx_mdl},'_largeScale_noLocs');
        if(do_save_fig)
            export_fig (fullfile(save_fig_path,figName),'-painters','-transparent','-pdf')
        end
    end
end

