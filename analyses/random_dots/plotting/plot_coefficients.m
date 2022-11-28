function [ah] = plot_coefficients(ah, mdl, plot_pars)
%{ Helper function to plot the coefficients of the linear models 
%  that regress overlap between the resdyn and task planes, veruses
%  properties of the resdyn (Ev/rot_freq/sv). 
% Input
% - ah - axis handle
% - mdl (cell array) - of size n_align x n_response_vars. Each cell
%   contains the linear model fit of a specific response variable
%   (angle_with_x, where x = {choice,time,JPC12,jPC34}) agaisn the
%   predictor variables (typically properties of residual dynamics)
% - plot_pars (struct) - cntaing plotting params 
%
% Output
% -ah - axis handle
%}
ylabels = {};
for ivar = 1: length(plot_pars.coef_plot.coefs_to_plot)
    c_count = 1;
    for iev = 1: size(mdl,1)
        for i_rv = 1: size(mdl,2)
            
            mdl_to_plot =  mdl{iev, i_rv};
            [coef_cis] = coefCI(mdl_to_plot);
            coef_names = mdl_to_plot.CoefficientNames;
            if(any(ismember(coef_names, plot_pars.coef_plot.coefs_to_plot{ivar})))
                coef_ = mdl_to_plot.Coefficients.Estimate(strcmp(coef_names,...
                    plot_pars.coef_plot.coefs_to_plot{ivar}));
                coef_cis = coef_cis(strcmp(coef_names,plot_pars.coef_plot.coefs_to_plot{ivar}),:);
                hh = ploterr(coef_,c_count,num2cell(coef_cis,1),[]);
                for ihh = 1:length(hh)
                    set(hh(ihh),'Color',plot_pars.coef_plot.coef_cols(ivar,:))
                end            

                plot(ah,coef_, c_count ,'o','MarkerEdgeColor',plot_pars.coef_plot.coef_cols(ivar,:),'MarkerFaceColor',[ 1 1 1]);
            end
         
            ylabels{c_count} = strcat(plot_pars.alignment_labels{iev} , ' ,' , plot_pars.plane_names{i_rv});
            c_count = c_count + 1;    
            
        end
    end
end
plot([0 0], [0.5 c_count - 1 + 0.5], '--k')
axis tight
set(ah,'ylim',[0.5 c_count - 1 + 0.5],'ytick',1:8,'yticklabels',ylabels);
xlabel(ah,'coeffcient magnitude')
ah = gca;
end