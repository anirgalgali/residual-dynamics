function [H] = plotXvTime(X,time,pars,varargin)

% X can either be a 2d array or 3d array (tensor) with the second dimension 
% having a size equal to the number of time-points to plot in the time series.
% The first dimension of X gets plotted according to different lineColors
% specified in pars. The third dimensions of X gets plotted using different
% line styles specified in pars

if(nargin > 3)
   h = varargin{1};
   ah = varargin{2};
%    hold(ah);
elseif(nargin == 3)
    
    figure; ah = gca; hold(ah,'plotbox',[1 1 1]);
    h = gcf;
end

set(h,'Position',pars.figPosition);

% if(size(pars.Color,1) == 1)
% 
%     line_cols = repmat(pars.Color,[size(X,1) 1]);
%     
% elseif(size(pars.Color,1) == size(X,1))
%     
%     line_cols = pars.Color;
%   
% end



for jj = 1:size(X,3)
    
    h_l = plot(ah,time,squeeze(X(:,:,jj))',pars.lineStyle{jj},'linewidth',pars.linewidth);
    set(h_l, {'color'}, num2cell(pars.Color{jj},2))
end

if(isfield(pars,'xLim') && isfield(pars,'yLim'))
    
    set(ah,'xlim',pars.xLim,'ylim',pars.yLim);
    
elseif(isfield(pars,'xLim') && ~isfield(pars,'yLim'))
    
    set(ah,'ylim',[min(X(:)) max(X(:))]); 
    
else
    
    axis tight;

end
set(ah,'plotbox',[1 1 1]);

H.fig_handle = h;
H.axes_handle = ah;

end