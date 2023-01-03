function [h,ah] = plotFlowField_v2(locs,flow,pars,varargin)
%{ This function plots the phase flow (flow) at locations on a 2D grid
%  (specified in locs). 
%  Inputs
% - locs (3D tensor) - the second and third dimension index the grid
%   locations. The first dimension specfies the (x,y) of each grid loc.
% - flow (3D tensor) - The corresponding flow vector at locations specified
%   in locs. The first dimension specifies the x and y components of the
%   flow.
% - pars (struct) - specifies the properties of the plots
% - varargin
%   - h : figure handle
%}  - ah : axes handle
if(nargin > 3)
    
   h = varargin{1};
   ah = varargin{2};hold all;
   
elseif(nargin == 3)
    
    figure; ah = gca; hold(ah,'plotbox',[1 1 1]);
    h = gcf;

end

if(isfield(pars,'figPosition'))
    set(h,'Position',pars.figPosition);
end

locs = locs(:,pars.start_X:pars.step_X:pars.end_X,pars.start_Y:pars.step_Y:pars.end_Y);
flow = flow(:,pars.start_X:pars.step_X:pars.end_X,pars.start_Y:pars.step_Y:pars.end_Y);

for jj = 1:size(flow,2)
    for kk = 1:size(flow,3)
        

        [arrow_shape_scaled] = createArrowParameters(pars.arrowShape,norm(pars.LineLength.*flow(:,jj,kk)));

        hh = arrows(locs(1,jj,kk),locs(2,jj,kk),...
            pars.LineLength.*flow(1,jj,kk),...
            pars.LineLength.*flow(2,jj,kk),...
            arrow_shape_scaled,'Cartesian','LineWidth',pars.lineWidth,'EdgeColor','none');
        
        set(hh,'Parent',ah);
        
        if(pars.plot_locs)
            plot(squeeze(locs(1,jj,kk)),squeeze(locs(2,jj,kk)),'o','Color',pars.dotColor,'markerfacecolor',pars.dotColor,'markersize',pars.dotSize);

        end
    end
end
%%
if(isfield(pars,'xLim'))
    set(ah,'xlim',pars.xLim)
end
if(isfield(pars,'yLim'))
    set(ah,'ylim',pars.yLim)
end

if(isfield(pars,'yTick'))
    set(ah,'ytick',pars.yTick)
end

if(isfield(pars,'xTick'))
    set(ah,'xtick',pars.xTick)
end

if(isfield(pars,'xTickLabel'))
    set(ah,'xticklabel',pars.xTickLabel)
end

if(isfield(pars,'yTickLabel'))
    set(ah,'yticklabel',pars.yTickLabel)
end


% flow_norms = sqrt(sum(abs(flow(:,:)).^2,1));
% ht = text(pars.xLim(2),pars.yLim(2),sprintf('%.3f',pars.LineLength*max(flow_norms)));
% set(ht,'horizontalalignment','right','verticalalignment','bottom');
end