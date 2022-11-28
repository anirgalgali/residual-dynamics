function [data_out,tLims_ind] = extractAnalysisTimeWindow(data,analysis_window)
%{ Extracts a time-window of the data for further analyses
% Input
% - data(struct) - containing responses and associated task/trial variables
% - analysis_window(struct) - containing two fields - tStart and tEnd,
% indicating the start and end times of the analysis window
% Output
% - data_out (struct) - filtered data belonging to the specified window
% - tLims_ind (array, logical) - of size n_times x 1, 1s indicate time
%   points that are selected
% Author : Aniruddh Galgali (Jan 2018). . Adapted from the TDR toolbox
% developed by Valerio Mante .
%}

data_out = data;
analysis_window_start = analysis_window.tStart;
analysis_window_end = analysis_window.tEnd;

tLims_ind = data.time >= analysis_window_start & data.time <= analysis_window_end;


data_out.time = data.time(tLims_ind);
data_out.response = data.response(:,tLims_ind,:);
if(isfield(data_out,'task_event'))
    data_out.task_event_align.time_window(1) = analysis_window_start;
    data_out.task_event_align.time_window(2) = analysis_window_end;
end

if(isfield(data_out,'time_rel'))
    data_out.time_rel = data.time_rel(tLims_ind);
end

if(isfield(data_out,'time_iev'))
    data_out.time_iev = data.time_iev(tLims_ind);
    
end

end