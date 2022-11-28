function [ah] = formatXvsTimeplot(ah,timeOfEvents,plot_stability,varargin)
%{ Helper function to format the plot of eige/signular values versus time.
% Inputs
% -ah - axis handle
% - timeOfEvents(array) - indicating the time-bins that correspond to event
%    onset/offset
% - plot_stability (boolean) - to indicate EV = 1 on the plots
%}
yLims = get(ah,'ylim');
xLims = get(ah,'xlim');
if(plot_stability)
    plot([xLims(1) xLims(2)],[1 1],'--','Color',[0.5 0.5 0.5]);
end

plotEventMarkers(ah,timeOfEvents,'x')

if(nargin == 4)

    snapshot_ticks = varargin{1};
    
    for hh = 1: length(snapshot_ticks)
        
        plot([snapshot_ticks(hh) snapshot_ticks(hh)],...
            [yLims(2) yLims(2) - 0.05],'-','Color',[0.5 0.5 0.5])
        
    end
    
end
 
ah = gca;

end

