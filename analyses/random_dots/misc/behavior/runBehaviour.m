%%
clearvars -except DIRS
clc
%%

doExtract = true;
if(doExtract)
    doSegment = true;
    extractBehaviour('Vito',doSegment);
end
    %%
close all
clear plotpars
line_col = cbrewer('qual','Dark2',8);
plotpars.marker_cols = repmat(line_col(1,:),[4 1]);
plotpars.marker_type = {'^','d','o','s'};
plotpars.marker_line_type = {'-^','-d','-o','-s'};
plotpars.markersize = 20;
plotpars.marker_face_alpha = 0.4;
plotpars.marker_edge_alpha = 0;
plotpars.coh_bins = [-80 -30 -10 0 10 30 80]; 
plotpars.line_type = 'fit_logistic';
plotpars.markercol = plotpars.marker_cols;
plotpars.plotAllData = true;
plotpars.plotBinnedAverage = false;
plotpars.fp_markersize = 8;
plotpars.fp_col = plotpars.marker_cols;
plotpars.targ_markersize = 3;
plotpars.pref_col = [0 0 1];
plotpars.anti_col = [1 0 0];
plotpars.FontSize = 8;
plotpars.FontName = 'Helvetica';
plotpars.save_fig_path = './figures/'; % enter a valid save path
plotpars.doSave = false;
animals = {'Tex';'Vito'};
data_path = './data/analyses/';

if(isempty(strfind(DIRS.analysis,'neuraldata')))
    DIRS.analysis = fullfile(DIRS.analysis,'/neuraldata/');
end

for iani = 1:length(animals)
    plotBehaviour(animals{iani}, data_path, plotpars);
end
%%

