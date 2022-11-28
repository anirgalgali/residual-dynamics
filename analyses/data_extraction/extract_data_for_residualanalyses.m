%{ This script extracts the raw data and performs a bunch of pre-processing 
%  steps to prepare the data for the residual analyses.
%  1) Removes silent and non-stationary units
%  2) Segments sessions into smaller experiments
%  3) Bins responses at a specific temporal resolution.
%
%  Author: Aniruddh Galgali
%}
clearvars -except DIRS
clc
rseed = rng('default');
close all

%%

pexp_data = ProjectLoad('array_TDR_dotsTask_data');
animal = 'Tex';

switch animal
    case 'Tex'        
        sessions_to_exclude= {'2008_07_22T'};
    case 'Vito'
        sessions_to_exclude = {'2007_08_03V';'2007_01_19V'};
end


pars.root_project = pexp_data;
pars.session.exclude_list = sessions_to_exclude;
pars.session.target_discrepancy_threshold = 0.4;
pars.session.segmentation.nChangePts = 5;
pars.session.segmentation.changeCountThreshold = 10;
pars.session.segmentation.thresholdRelMeanChange = 0.15;
pars.session.segmentation.numTrialAvgWin = 50;

pars.align.Fs = 200;
pars.align.aligns_to_keep = {'stimoncd';'outfixcd'};
pars.align.pre_analysis_window.tLims{1} = [-0.5 1.5];
pars.align.pre_analysis_window.tLims{2} = [-1.2 0.8];

pars.unit_exclusion.do_remove_quiescent = true;
pars.unit_exclusion.do_remove_lowchoiceselective = false;
pars.unit_exclusion.quiescentUnits.threshold = 1;
pars.unit_exclusion.choiceSel.threshold = 2;
pars.unit_exclusion.choiceSel.numTimeSegments = 5;
pars.unit_exclusion.do_remove_nonstat = true;
pars.unit_exclusion.nonStat.winStep = 20;
pars.unit_exclusion.nonStat.winLength = 50;
pars.unit_exclusion.nonStat.minTrialsToAvg = 10;
pars.unit_exclusion.nonStat.thresholdToUse = 'mean_std'; % always this option
pars.unit_exclusion.nonStat.minTrialsinSeg = 200;
pars.unit_exclusion.nonStat.percentileThreshold = 99;
pars.unit_exclusion.nonStat.nTrialsLookAhead = 50;
pars.unit_exclusion.nonStat.relChangeThreshold = 0.2;
assert(length(pars.align.pre_analysis_window.tLims) == length(pars.align.aligns_to_keep))

pars.single_trial.analysis_window.tLims{1} = [-0.2 1.0];
pars.single_trial.analysis_window.tLims{2} = [-0.7 0.5];
pars.single_trial.bin_size = 90;
pars.single_trial.use_sqrt_transform = true;
pars.single_trial.use_rates = false;
pars.single_trial.delay_length_bin_boundaries = [0 0.4 0.6 0.8 1.0 1.5];

project_name = sprintf('%s%d%s','array_reppasdotsTask_datasegmented_binsize=', pars.single_trial.bin_size,'ms');
pexp_segmented = ProjectCreate(project_name);

save_project.name = project_name;
save_project.root_dir = DIRS.analysis;
save_project.loader = pexp_segmented;

extract_and_preprocess_data(pexp_data, animal, save_project, pars);
