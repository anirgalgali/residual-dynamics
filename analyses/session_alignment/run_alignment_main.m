%{ This script performs the session alignment and computes the residuals.
%  and saves the results to disk.
%  Author: Aniruddh Galgali
% Note: You may want to change the paths to point to where the data is stored.
%}
%%
clearvars -except DIRS
clc
%% Choosing animal

animal = 'Tex';

%% Extracting alignment and residuals for a bin-size of 45ms
project_name = 'array_reppasdotsTask_datasegmented_binsize=45ms';
pexp_data = ProjectLoad_v2(project_name);
clear pars_
pars_.alignment.nPCs_align = 20;
pars_.alignment.type = 'full';

if(strcmp(pars_.alignment.type,'precomputed'))
   file_path = '';
   file_name = '';
end

pars_.align.condition_vars = {'targ_dir'};
pars_.align_proj.condition_vars = {'targ_dir'};

pars_.residual.condition_vars = {{'targ_dir','dots_coh'},{'targ_dir','delay_length'}};
pars_.residual.exclude_trials = {{'delay_length', 1}};

pars_.config_vars = 'target_config';
[data_aligned] =  compute_alignment_and_residuals(pexp_data, animal, pars_);
cd('/Users/Aniruddh/Work_PhD/Residual Dynamics/residual_dynamics');
save_file_path = './data/processedneuraldata/';
save_file_name = sprintf('%s%s',animal,'_aligned_reppasdotsTask_binsize=45ms.mat');
save(fullfile(save_file_path,save_file_name),'data_aligned','-v7.3');
%% Extractiing non-aligned residuals (for the single session fits)

project_name = 'array_reppasdotsTask_datasegmented_binsize=45ms';
pexp_data = ProjectLoad_v2(project_name);
clear pars_
pars_.alignment.nPCs_align = 20;
pars_.alignment.type = 'none';

if(strcmp(pars_.alignment.type,'precomputed'))
   file_path = '';
   file_name = '';
end

pars_.align.condition_vars = {'targ_dir'};
pars_.align_proj.condition_vars = {'targ_dir'};

pars_.residual.condition_vars = {{'targ_dir','dots_coh'},{'targ_dir','delay_length'}};
pars_.residual.exclude_trials = {{'delay_length', 1}};

pars_.config_vars = 'target_config';
[data_aligned] =  compute_alignment_and_residuals(pexp_data, animal, pars_);
save_file_path = '/Volumes/WD Passport/WorkPC_2021Backup_final/Projects/ResDynPFC/analysis_datasets/neuraldata/';
save_file_name = sprintf('%s%s',animal,'_nonaligned_reppasdotsTask_binsize=45ms.mat');
save(fullfile(save_file_path,save_file_name),'data_aligned','-v7.3');
%% Extracting alignment and residuals for a bin-size of 15ms

% This uses the precomputed alignment for the 45ms bin-size
clear pars_
pars_.alignment.nPCs_align = 20;
pars_.alignment.type = 'precomputed';
project_name = 'array_reppasdotsTask_datasegmented_binsize=15ms';
pexp_data = ProjectLoad_v2(project_name);

if(strcmp(pars_.alignment.type,'precomputed'))
   cd('/Users/Aniruddh/Work_PhD/Residual Dynamics/residual_dynamics');
   file_path = './data/processedneuraldata/';
   file_name = sprintf('%s%s',animal,'_aligned_reppasdotsTask_binsize=45ms.mat');
   load(fullfile(file_path,file_name))
   pars_.precomputed_alignment.U_orth = cellfun(@(x) x.U_align,data_aligned.aligned_single_trials,'uni',false);
   pars_.precomputed_alignment.X_mu = cellfun(@(x) x.X_mu, data_aligned.aligned_single_trials,'uni',false);
end

pars_.align.condition_vars = {'targ_dir'};
pars_.align_proj.condition_vars = {'targ_dir'};

pars_.residual.condition_vars = {{'targ_dir','dots_coh'},{'targ_dir','delay_length'}};
pars_.residual.exclude_trials = {{'delay_length', 1}};

pars_.config_vars = 'target_config';
[data_aligned] =  compute_alignment_and_residuals(pexp_data, animal, pars_);
save_file_path = '/Volumes/WD Passport/WorkPC_2021Backup_final/Projects/ResDynPFC/analysis_datasets/neuraldata/';
if(strcmp(pars_.alignment.type,'precomputed'))
    save_file_name = sprintf('%s%s',animal,'_aligned_reppasdotsTask_binsize=15ms_alignbw=45ms.mat');
elseif(strcmp(pars_.alignment.type,'full'))
    save_file_name = sprintf('%s%s',animal,'_aligned_reppasdotsTask_binsize=15ms.mat');   
end
save(fullfile(save_file_path,save_file_name),'data_aligned','-v7.3');
%% Extracting alignment and residuals for a bin-size of 30ms

% This uses the precomputed alignment for the 45ms bin-size
clear pars_
pars_.alignment.nPCs_align = 20;
pars_.alignment.type = 'precomputed';
project_name = 'array_reppasdotsTask_datasegmented_binsize=30ms';
pexp_data = ProjectLoad_v2(project_name);

if(strcmp(pars_.alignment.type,'precomputed'))
   cd('/Users/Aniruddh/Work_PhD/Residual Dynamics/residual_dynamics');
   file_path = './data/processedneuraldata/';
   file_name = sprintf('%s%s',animal,'_aligned_reppasdotsTask_binsize=45ms.mat');
   load(fullfile(file_path,file_name))
   pars_.precomputed_alignment.U_orth = cellfun(@(x) x.U_align,data_aligned.aligned_single_trials,'uni',false);
   pars_.precomputed_alignment.X_mu = cellfun(@(x) x.X_mu, data_aligned.aligned_single_trials,'uni',false);
end

pars_.align.condition_vars = {'targ_dir'};
pars_.align_proj.condition_vars = {'targ_dir'};

pars_.residual.condition_vars = {{'targ_dir','dots_coh'},{'targ_dir','delay_length'}};
pars_.residual.exclude_trials = {{'delay_length', 1}};

pars_.config_vars = 'target_config';
[data_aligned] =  compute_alignment_and_residuals(pexp_data, animal, pars_);
save_file_path = '/Volumes/WD Passport/WorkPC_2021Backup_final/Projects/ResDynPFC/analysis_datasets/neuraldata/';
if(strcmp(pars_.alignment.type,'precomputed'))
    save_file_name = sprintf('%s%s',animal,'_aligned_reppasdotsTask_binsize=30ms_alignbw=45ms.mat');
elseif(strcmp(pars_.alignment.type,'full'))
    save_file_name = sprintf('%s%s',animal,'_aligned_reppasdotsTask_binsize=30ms.mat');   
end
save(fullfile(save_file_path,save_file_name),'data_aligned','-v7.3');
%% Extracting alignment and residuals for a bin-size of 60ms

% This uses the precomputed alignment for the 45ms bin-size
clear pars_
pars_.alignment.nPCs_align = 20;
pars_.alignment.type = 'precomputed';
project_name = 'array_reppasdotsTask_datasegmented_binsize=60ms';
pexp_data = ProjectLoad_v2(project_name);

if(strcmp(pars_.alignment.type,'precomputed'))
   cd('/Users/Aniruddh/Work_PhD/Residual Dynamics/residual_dynamics');
   file_path = './data/processedneuraldata/';
   file_name = sprintf('%s%s',animal,'_aligned_reppasdotsTask_binsize=45ms.mat');
   load(fullfile(file_path,file_name))
   pars_.precomputed_alignment.U_orth = cellfun(@(x) x.U_align,data_aligned.aligned_single_trials,'uni',false);
   pars_.precomputed_alignment.X_mu = cellfun(@(x) x.X_mu, data_aligned.aligned_single_trials,'uni',false);
end

pars_.align.condition_vars = {'targ_dir'};
pars_.align_proj.condition_vars = {'targ_dir'};

pars_.residual.condition_vars = {{'targ_dir','dots_coh'},{'targ_dir','delay_length'}};
pars_.residual.exclude_trials = {{'delay_length', 1}};

pars_.config_vars = 'target_config';
[data_aligned] =  compute_alignment_and_residuals(pexp_data, animal, pars_);
save_file_path = '/Volumes/WD Passport/WorkPC_2021Backup_final/Projects/ResDynPFC/analysis_datasets/neuraldata/';
if(strcmp(pars_.alignment.type,'precomputed'))
    save_file_name = sprintf('%s%s',animal,'_aligned_reppasdotsTask_binsize=60ms_alignbw=45ms.mat');
elseif(strcmp(pars_.alignment.type,'full'))
    save_file_name = sprintf('%s%s',animal,'_aligned_reppasdotsTask_binsize=60ms.mat'); 
end
save(fullfile(save_file_path,save_file_name),'data_aligned','-v7.3');
%% Extracting alignment and residuals for a bin-size of 90ms

% This uses the precomputed alignment for the 45ms bin-size
clear pars_
pars_.alignment.nPCs_align = 20;
pars_.alignment.type = 'precomputed';
project_name = 'array_reppasdotsTask_datasegmented_binsize=90ms';
pexp_data = ProjectLoad_v2(project_name);

if(strcmp(pars_.alignment.type,'precomputed'))
   cd('/Users/Aniruddh/Work_PhD/Residual Dynamics/residual_dynamics');
   file_path = './data/processedneuraldata/';
   file_name = sprintf('%s%s',animal,'_aligned_reppasdotsTask_binsize=45ms.mat');
   load(fullfile(file_path,file_name))
   pars_.precomputed_alignment.U_orth = cellfun(@(x) x.U_align,data_aligned.aligned_single_trials,'uni',false);
   pars_.precomputed_alignment.X_mu = cellfun(@(x) x.X_mu, data_aligned.aligned_single_trials,'uni',false);
end

pars_.align.condition_vars = {'targ_dir'};
pars_.align_proj.condition_vars = {'targ_dir'};

pars_.residual.condition_vars = {{'targ_dir','dots_coh'},{'targ_dir','delay_length'}};
pars_.residual.exclude_trials = {{'delay_length', 1}};

pars_.config_vars = 'target_config';
[data_aligned] =  compute_alignment_and_residuals(pexp_data, animal, pars_);
save_file_path = '/Volumes/WD Passport/WorkPC_2021Backup_final/Projects/ResDynPFC/analysis_datasets/neuraldata/';
if(strcmp(pars_.alignment.type,'precomputed'))
    save_file_name = sprintf('%s%s',animal,'_aligned_reppasdotsTask_binsize=90ms_alignbw=45ms.mat');
elseif(strcmp(pars_.alignment.type,'full'))
    save_file_name = sprintf('%s%s',animal,'_aligned_reppasdotsTask_binsize=90ms.mat');  
end
save(fullfile(save_file_path,save_file_name),'data_aligned','-v7.3');