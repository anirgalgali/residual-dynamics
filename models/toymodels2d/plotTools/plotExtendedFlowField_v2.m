function [H] = plotExtendedFlowField_v2(U,V,pars,varargin)
    
if(nargin > 3)
   h = varargin{1};
   ah = varargin{2};hold all;
   
elseif(nargin == 3)
    figure; ah = gca; hold(ah,'plotbox',[1 1 1]);
    h = gcf;
end

set(h,'Position',pars.figPosition);



U = U(:,pars.start_X:pars.step_X:pars.end_X,pars.start_Y:pars.step_Y:pars.end_Y);
V = V(:,pars.start_X:pars.step_X:pars.end_X,pars.start_Y:pars.step_Y:pars.end_Y);

for jj = 1:size(V,2)
    for kk = 1:size(V,3)
        


        plot(ah,squeeze(U(:,jj,kk)),squeeze(V(:,jj,kk)),'Color',pars.Color)
        
%         hh = arrows(U(pars.arrow_start_idx,jj,kk),V(pars.arrow_start_idx,jj,kk),...
%             pars.LineLength.*(U(pars.arrow_start_idx + 1,jj,kk) - U(pars.arrow_start_idx,jj,kk)),...
%             pars.LineLength*(V(pars.arrow_start_idx + 1,jj,kk) - V(pars.arrow_start_idx,jj,kk)),...
%             pars.arrowShape,'Cartesian','LineWidth',pars.lineWidth);
%         
%         set(hh,'Parent',ah);
        

    end
end
%%
set(ah,'xlim',pars.xLim,'ylim',pars.yLim,'xtick',pars.xTick,'ytick',pars.yTick,'xticklabel',pars.xTickLabel,'yticklabel',pars.yTickLabel);    
H.fig_handle = h;
H.axes_handle = ah;
end