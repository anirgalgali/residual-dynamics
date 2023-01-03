function [A_d] =  discretize_dynamics_matrix(A, dt)
%{
This function discretizes the dynamics matrix that specifies a continuous
time state-space model.

Inputs : A - matrix/tensor
            (latent_dim x latent_dim) or (latent_dim x latent_dim x num_times)
            A matrix or tensor specifying the time-invariant or time-varying
            dynamics matrix. 
         dt - scalar
            step-size for discretization

Outputs : A_d - matrix/tensor
            (latent_dim x latent_dim) or (latent_dim x latent_dim x num_times)
            A matrix or tensor specifying the "discretized" time-invariant or
            time-varying dynamics matrix. 

%}
if(ismatrix(A))
    
   A_d = expm(A.*dt);
    
else
    n_times = size(A,3);
    A_d = NaN(size(A));
    for tt = 1:n_times

        A_d(:,:,tt) = expm(A(:,:,tt).*dt);

    end
    
end

end