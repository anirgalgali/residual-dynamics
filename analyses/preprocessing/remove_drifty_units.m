function [drifty_idxs,fr_stats] = remove_drifty_units(spike_fr,pars)

% spike_fr = data.spike_rate_temporalAverage;
% removed_indices = data.removed_units;
[fr_stats] = fit_linear_firing_rate(spike_fr);
stat_names = fieldnames(fr_stats);
assert(isfield(fr_stats, pars.stat_type),'wrong statistic');

stat = fr_stats.(pars.stat_type);
drifty_idxs = stat >= pars.thresh_stat;

[~,idx] = sort(stat,'descend');
for i_stat = 1: length(stat_names)
    fr_stats.(stat_names{i_stat}) = fr_stats.(stat_names{i_stat})(idx);
end

fr_stats.resorted_unit_idxs = idx;



end

function [fr_stats] = fit_linear_firing_rate(y)

[num_units,num_trials] = size(y);
r2 = zeros(num_units,1);
slopes = zeros(num_units,1);
slope_95CI = zeros(num_units,1);
intercepts = zeros(num_units,1);
F_stat = zeros(num_units,1);
pVal = zeros(num_units,1);
X = [1: num_trials]';

for i_unit = 1: num_units
    
    mdl = fitlm(X,y(i_unit,:)','robustOpts','huber');
    r2(i_unit) = mdl.Rsquared.Ordinary;
    slopes(i_unit) = mdl.Coefficients.Estimate(2);
    slope_95CI(i_unit) = 1.96*(mdl.Coefficients.SE(2));
    intercepts(i_unit) = mdl.Coefficients.Estimate(1);
    pVal(i_unit) = mdl.Coefficients.pValue(2);
    mdl_an = anova(mdl);
    F_stat(i_unit) = mdl_an.F(1);
    
end

fr_stats.rsquared = r2;
fr_stats.coeff = slopes;
fr_stats.intercept = intercepts;
fr_stats.coeff_CI = slope_95CI;
fr_stats.F_val = F_stat;
fr_stats.coeff_pVal = pVal;


end