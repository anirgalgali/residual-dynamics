function [tmax] = extract_timeidx_of_maxev(result)
%{ This function extracts the time-index corresponding to the largest
%  eigenvalue over time
%  INPUT
%  - data : struct containing the results of the residual dynamics fit
%  OUTPUT
%  - tmax : (n_conds x 1) - The tmax corresponding to each task condition
%  (choices) for which residual dynamics is fitted
%  Author: Aniruddh Galgali (Nov 2020)
%}
n_choices = length(result.final_model);
tmax = zeros(length(result.final_model),1);
for icond = 1: n_choices
    
    max_ev_t = max(result.final_model{icond}.eigVal,[],1);
    max_ev_delta = max(max_ev_t) - min(max_ev_t);
    
    % This value of 0.05 establishes the threshold of change in the largest
    % and smallest value of the maximum ev across time that is considered 
    % to be signficant. If the difference is lesse than 0.05, then the time
    % course of the maximum ev is just considered to be time-invariant, and
    % the last time index is returned
    
    if(max_ev_delta <= 0.05) 
         tmax(icond) = length(max_ev_t);
    else

        [~,tmax(icond)] = max(max_ev_t);

    end
    
end
end