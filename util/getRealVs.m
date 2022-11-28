function Vr = getRealVs(V,evals)

% Function to extract a pair of real eigenvectors from complex eigenvectors
% Input 
% - V (matrix) - of size n_dim x n_dim, each column corresponds to a complex
%   eigenvector.
% - evals (array) - of size n_dim x 1, each value is the corresponding
%   complex-valued eigenvalue.
% Outputs
% - Vr (matrix) - of size n_dim x n_dim, each column corresponds to a real
%   eigenvector.
%
% Author : Mark Churchland/John Cunninghma. Directly adapted from the jPCA toolbox
% get real vectors made from the eigenvectors

% by paying attention to this order, things will always rotate CCW
if imag(evals(1))>0  % if the eigenvalue with negative imaginary component comes first
    Vr = [V(:,1) + V(:,2), (V(:,1) - V(:,2))*1i];
else
    Vr = [V(:,2) + V(:,1), (V(:,2) - V(:,1))*1i];
end
if(any(imag(Vr) ~= 0))
    Vr = normc(real(Vr));
else
    Vr = normc(Vr);
end

end