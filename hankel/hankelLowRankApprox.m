function [Hr] = hankelLowRankApprox(H,r)
%{ 

Computes a (hard thresholded) r-rank approximation of the hankel matrix H.
H can either be asingle hankel matrix or a temporal sequence of hankel matrices.

Inputs : H - either a 2D matrix  of size n_dim_hankel x n_dim_hankel,  or a
             3D tensor of size n_dim_hankel x n_dim_hankel x n_times_hankel
         r - scalar specifying the rank 

Outputs: Hr - same size as H, but instead each matrix is of rank-r

        Author : Aniruddh Galgali, Dec 2017

%}
if(ndims(H) == 3)
    Hr = NaN(size(H));
    for tt = 1:size(H,3)
        
        [U,S,V] = svd(H(:,:,tt));
        Hr(:,:,tt) = U(:,1:r)*S(1:r,1:r)*V(:,1:r)';
        
    end
    
elseif(ndims(H) == 2)
    
    [U,S,V] = svd(H);
    Hr = U(:,1:r)*S(1:r,1:r)*V(:,1:r)';
    
end

end