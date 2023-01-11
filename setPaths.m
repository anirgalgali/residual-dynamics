% add paths to source and auxilliary code toolboxes

% This should be the path to where you downloaded the Github repo
addpath(genpath('/Users/Aniruddh/Work_PhD/Residual Dynamics/residual_dynamics/'))

% set paths to neural data 
global DIRS

% this is the path to the raw neural data (this may vary depending on where
% one stores the raw neural data
if isunix
    % This is the path to the folder where the data has been stored.
    % Typically, the data is stored in a folder titled 'data' within the
    % main code directory. Change this accordingly
    DIRS.analysis = '/Volumes/WD Passport/data';
else
    DIRS.analysis = '';
end
DIRS.user = 'AG';

% set plotting defaults
set(groot, ...
'DefaultFigureColor', 'w', ...
'DefaultAxesLineWidth', 0.5, ...
'DefaultAxesXColor', 'k', ...
'DefaultAxesYColor', 'k', ...
'DefaultAxesFontUnits', 'points', ...
'DefaultAxesFontSize', 8, ...
'DefaultAxesFontName', 'Helvetica', ...
'DefaultLineLineWidth', 1, ...
'DefaultTextFontUnits', 'Points', ...
'DefaultTextFontSize', 12, ...
'DefaultTextFontName', 'Helvetica', ...
'DefaultAxesBox', 'off', ...
'DefaultAxesTickLength', [0.02 0.025]);
 
% set the tickdirs to go out - need this specific order
set(groot, 'DefaultAxesTickDir', 'out');
set(groot, 'DefaultAxesTickDirMode', 'manual');

% Change to current working directory

% This should be the path to where you downloaded the Github repo
cd('/Users/Aniruddh/Work_PhD/Residual Dynamics/residual_dynamics/')