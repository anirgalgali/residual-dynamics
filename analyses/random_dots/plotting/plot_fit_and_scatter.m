function [ah] = plot_fit_and_scatter(ah, x, y, xvar_label, mdls, plotpars, varargin)
%{ Helper function that plots the scatter and corresponding linear fits 
% that show the relationship between the overlap of the resdyn and task 
% and properties of the resdyn (Ev/rot_freq/sv). 
% Input
% - ah - axis handle
% - x(array) - values of predcitor variable to plot
% - y(array) - values of response variable to plots
% - xvar_label - string indicating the identity of the predcitor variable
% - mdls (native lm model class) - the linear model obtained after
%   regressing y v x.
% - plot_pars (struct) - cntaing plotting params 
% - varargin
%   - mean(struct) - means of the redictor variables (used in case
%     predictors were centred before performing linear regression)
%
% Output
% -ah - axis handle
%}  

if(nargin > 6) 
   xmean = varargin{1}; 
else
    xmean = 0;
end

plot(ah, x + xmean, y,'.','markersize',plotpars.scatter.marker_size,'color',plotpars.scatter.marker_col)
set(ah,'ylim',[0 90],'ytick',[0:30:90],'yticklabel',[0:30:90])

xx = linspace(min(x), max(x), 100);


coef_names = mdls.CoefficientNames;

if(any(ismember(coef_names, xvar_label)))
    coef1_ = mdls.Coefficients.Estimate(strcmp(coef_names, xvar_label));
else
    coef1_ = 0.0;
end

coef0_ = mdls.Coefficients.Estimate(strcmp(coef_names,'(Intercept)'));
plot(xx + xmean, (coef0_ - coef1_.*xmean) + (coef1_.* (xx + xmean)), plotpars.scatter.mdl_fit_style, 'Color',[0 0 0])
  
xlabel(ah,sprintf('%s%s%s','|',xvar_label,'|'))
   
end