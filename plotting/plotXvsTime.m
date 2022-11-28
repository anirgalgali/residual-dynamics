function [varargout] = plotXvsTime(ah,data,time,plotpars,varargin)
%{ Function thata makes a pretty plot of quantity in 'data' versus 'time'
% Inputs
% - ah - axis handle
% - data(cell array) - each element of the cell aray corresponds to a lline
%   in the plot.
% - time (array) - time bins
% - plotpars (struct) - high level plotting parameters.
%}
% data should be a cell array. each cell should contain a time-series you
% want to plot within the same figure

% Setting default plot parameters


time_rel  = time;
time_iev = ones(length(time),1);
ci = [];
assignopts(who, varargin);

hold all;
if(isfield(plotpars,'FontSize'))
    font_size = plotpars.FontSize;
else
    font_size = 8;
end

if(isfield(plotpars,'FontName'))
    font_name = plotpars.FontName;
else
    font_name = 'Helvetica';
end

n_x = length(data);

unique_time_ievs = unique(time_iev);

for i_x = 1:n_x
    
    for kk = 1: length(unique_time_ievs)
        
        
        if(isempty(ci))
        
            h_l = plot(time(time_iev == unique_time_ievs(kk)), ...
                data{i_x}(time_iev == unique_time_ievs(kk)),plotpars.line_style{i_x});
            
            set(h_l,'Color',plotpars.col_lines(i_x,:),'linewidth',plotpars.linewidth);
            set(h_l,'Parent',ah);
        
        else
           
            switch plotpars.ci_type
                
                case 'patch'
            
                    yy = data{i_x}(time_iev == unique_time_ievs(kk));
                    yy_err = ci{i_x}(time_iev == unique_time_ievs(kk),:);
                    nan_idxs = isnan(yy_err(:,1));
                    time_ev = time(time_iev == unique_time_ievs(kk));

                    
                    [h_l,h_p] = boundedline(time_ev(~nan_idxs), ...
                        yy(~nan_idxs),yy_err(~nan_idxs,:),plotpars.line_style{i_x});
                    
                    if(plotpars.do_color_edge)
                        set(h_p,'FaceColor',plotpars.col_lines(i_x,:),...
                            'facealpha',plotpars.face_alpha, 'edgecolor',plotpars.col_lines(i_x,:));
                        
                    else
                        set(h_p,'FaceColor',plotpars.col_lines(i_x,:),...
                            'facealpha',plotpars.face_alpha, 'edgecolor','none');
                        
                    end
 
                    set(h_l,'Color',plotpars.col_lines(i_x,:),'linewidth',plotpars.linewidth);
                    set(h_l,'Parent',ah);
                    
                case 'bar'
                    
                    yy = data{i_x}(time_iev == unique_time_ievs(kk));
                    time_ev = time(time_iev == unique_time_ievs(kk));

                    
                    ci_lower = ci{i_x}(time_iev == unique_time_ievs(kk),1);
                    ci_upper = ci{i_x}(time_iev == unique_time_ievs(kk),2);

                    yy_err = {ci_lower,ci_upper};
                    h_l = ploterr(time_ev,yy,[],yy_err,'-','abshhy', plotpars.bar_length);
                    set(h_l(1),'Color',plotpars.col_lines(i_x,:),'linewidth',plotpars.linewidth);
                    set(h_l(2),'Color',plotpars.col_lines(i_x,:),'linewidth',plotpars.bar_width);
                    
            end
            
            set(h_l,'Parent',ah);
                
        end
        
               
        
        if(plotpars.do_plot_markers)
            
            h_m = plot(ah,time( time_iev == unique_time_ievs(kk)), ...
                data{i_x}(time_iev == unique_time_ievs(kk)),plotpars.marker_style(i_x,:));
            
            set(h_m,'markeredgecolor',plotpars.col_lines(i_x,:) ,'markersize',plotpars.markersize,'markerfacecolor',plotpars.markerfacecol(i_x,:));
        end
    
        
        
    end
end
    
axis tight
if(isfield(plotpars,'ylim'))
    set(gca,'ylim',plotpars.ylim,'ytick',plotpars.ytick,'yticklabel',plotpars.yticklabel)
end


if(isfield(plotpars,'xLabel'))
    set(get(gca,'XLabel'),'String',plotpars.xLabel,'FontName',font_name,'FontSize',font_size);
end

if(isfield(plotpars,'yLabel'))
    set(get(gca,'YLabel'),'String',plotpars.yLabel,'FontName',font_name,'FontSize',font_size);
end

varargout{1} = gca;

end