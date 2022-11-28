function plotEventMarkers(ah,timeOfEvents,type)
%{ Helper function to plot ticks at time-bins indicating event onset/offset
% Inputs
% -ah - axis handle
% - timeOfEvents(array) - indicating the time-bins that correspond to event
%    onset/offset
% - type (string) - can be either {'x','y' 'xy'}, indicates the axes on
%   which to show the ticks
%}
yLims = get(ah,'ylim');
xLims = get(ah,'xlim');

switch type
    
    case 'x'
        for jj = 1:length(timeOfEvents)
            
            
            plot(ah,[timeOfEvents(jj) timeOfEvents(jj)],[yLims(1) yLims(1) + 1/20*(yLims(2) - yLims(1))],'-k','linewidth',2);
            
        end
        
    case 'xy'
        
        for jj = 1:length(timeOfEvents)
            
            
            plot(ah,[timeOfEvents(jj) timeOfEvents(jj)],[yLims(1) yLims(1) + 1/20*(yLims(2) - yLims(1))],'-k','linewidth',2);
            plot(ah,[xLims(1) xLims(1) + 1/20*(yLims(2) - yLims(1))],[timeOfEvents(jj) timeOfEvents(jj)],'-k','linewidth',2);
        
        end
        
end
end