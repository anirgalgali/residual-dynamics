function [A] =  generateParametersforNonNormalDynamics(pars)

%{ This function returns the dynamics update equation corresponding to 
% a non-normal dynamics model
%}



% input_dir = [1 1; -1 1]'; % input directions
% input_sigma = [1 1]; % input noise standard deviation
% % input_sigma = [.3 .3]; % input noise standard deviation
% input_mean  = [0 0]; % input mean
% start_sigma = [.3 .3]; % starting point standard deviation
% start_mean  = [0 0]; % starting point mean



%--- Generate the interaction matrix ---
% 90 deg rotation matrix
Mrot = [0 -1; 1 0];

% Directions of first right and left eigenvectors


% Direction of second right eigenvector
r1 = Mrot*normc(pars.l0);

% Matrix or rights
R = [pars.r0 r1];

% Matrix of lefts
L = inv(R);

% Diagonal matrix of eigenvalues
E = diag(pars.dyn_eigenval);

% Interaction matrix
A = R*E*L;


end