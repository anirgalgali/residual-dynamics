function [timeOfEvents] = findEventMarkers(time_abs,time_rel,time_iev)
% Function to extract time bins that correspond to event onsets
% Input 
% - time_abs (array) - of size n_times x 1, indicates the absolute times
%   for all time-bins in the data (typically concatenated over all temporal
%   alignments)
% - time_rel (array) - of size n_times x 1, indicates the 'relative' times
%   in each temporal alignment.
% - time_labels (array) - of size n_times x 1, each element is a
%   numerical assignment of the time bin to a specific task epoch.
% Outputs
% - timeOfEvents (array) 
%
% Author : Aniruddh Galgali.

% saccade and stimulus onset
unique_time_iev = unique(time_iev);
% These are hard-coded. Although 0 is typically the onset of events, when
% aligning neural data to event onset. The value of 0.8 is specific to the
% dots task dataset, and may need to be changed for different datasets.
% TODO: Generalize this to arbitrary event onset/offset times.
time_ev = {[0;0.8];0}; 
ev_count = 1;
timeOfEvents = [];
for ii = 1:length(unique_time_iev)
    
    all_events = time_ev{ii};
    
    for jj = 1:length(all_events)
       
        timeOfEvents(ev_count) = interp1(time_rel(time_iev == unique_time_iev(ii)),...
            time_abs(time_iev == unique_time_iev(ii)),all_events(jj));
        
        ev_count = ev_count + 1;
        
    end
    
    
end



end