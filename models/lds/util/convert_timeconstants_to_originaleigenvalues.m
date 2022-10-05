function [x] = convert_timeconstants_to_originaleigenvalues(tau, dt)
% Converts time-constants into eigenvalues of the dynamics. 

% Inputs 

% tau - time-constant (positive)
% dt - time step in which the original parameters are specified

x = exp(-dt./tau);

end