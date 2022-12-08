function tdr()
% The TDR ('Targeted Dimensionality Reduction') toolbox implements the 
% dimensionality reduction methods described in the paper:
%
% 'Context-dependent computation by recurrrent dynamics in prefrontal cortex'
%  by Mante, Sussillo, Shenoy, and Newsome, Nature, 2013.
%
% The neurophysiology and behavioral datasets described in the paper are
% included in the toolbox.
%
% For questions email Valerio Mante (valerio@ini.phys.ethz.ch)
%
% The script 'pfc TDR analysis' loads the neurophysiology data, and calls
% the functions necessary to reproduce the main figures in the paper.
%
% Several TDR functions operate on data. The data comes in two formats:
%
% (1) Sequentially recorded responses. 
%     The raw data from the paper is saved in this format. Responses from 
%     the various units (single units and multi-units) are mostly recorded 
%     sequentially on separate days.
%
%     data.unit(i): the data for unit i. Has the following fields:
%        .response: binned spikes for unit i. An array of dimensions [ntr npt]
%            (ntr: number of trials, npt: number of time points)
%        .task_variable: a structure of task variables, each provided as a
%            vector of size [ntr 1]. For instance:
%        .task_variable.stim_dir contains the direction coherence of the
%            random dots on each trial.
%        .dimension: the unit label.
%     data.time: the time axis [1 npt], common to all units.
%
% (2) Simultanously recorded responses.
%     Data from simultaneously recorded units, or condition averaged
%     responses from either simultaously or sequentially recorded data.
%     
%     This data type is saved in a structure with the following fields:
%        .response: binned spikes for all units. An array of dimensions
%            [nun npt ncd] (nun: number of units, npt: number of time
%            points, ncd: number of trials/conditions).
%        .task_variable: a structure of task variables, each provided as a
%            vector of size [ncd 1]. For instance:
%        .task_variable.stim_dir contains the direction coherence of the
%            random dots on each trial. For instance, coherences of 
%            [-50 -20 -5 5 20 50 -50 50] become indeces [1 2 3 4 5 6 1 6].
%        .n_trial: for condition averaged responses, these are is the number 
%            of trials per conditions, a vector of size [1 ncd]. For
%            simultaneously recorded response this is row of ones.
%        .time: the time axis [1 npt].
%        .dimension: the labels of each dimension (e.g. 'unit_1'). Cell
%            array of dimension {nun 1}


































