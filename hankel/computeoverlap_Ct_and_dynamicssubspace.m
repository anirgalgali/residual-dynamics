function [overlap] =  computeoverlap_Ct_and_dynamicssubspace(C_t, U)
%{ 

Computes the subspace overlap between the final time-invariant dynamics
subspace (U) and the sequence of observation matrices (C_t) that went into
constructing the dynamics subspace

Inputs : U - dynamics subspace:  matrix of size n_dim x n_dim, where the 
             columns are ordered in decreasing order of importance.
         C_t - time-varying observation matrices : cell array of length
         (n_conds x n_hankel_times) x 1, where each entry corresponds the
         observation matrix at a particular time-point and  condition. 

Outputs: Hr - same size as H, but instead each matrix is of rank-r

        Author : Aniruddh Galgali, Dec 2021

%}

overlap = cell(size(U,2),1);
for idim = 1: size(U,2)
    
    overlap{idim} = cell2mat(cellfun(@(x) ...
        rad2deg(subspace(x, U(:,1:idim))),C_t,'uni',false));
end

end