function [tau] = convert_originaleigenvalues_to_timeconstants(x, dt)
% Converts the recovered parameters of the dynamics and noise into
% time-constants. 

%Use this function to obtain time constants from the recovered parameter
%vector from the fmincon optimization

% Inputs 

% x - recovered/optimized eigenvalues
% dt - time step in which the original parameters are specified

tau = abs(dt./log(x)); 

end