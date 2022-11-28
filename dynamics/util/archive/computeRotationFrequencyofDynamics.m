function [rot_freq] = computeRotationFrequencyofDynamics(A, dt, sort_type)

n_time = size(A,3);
n_dim = size(A,1);

switch sort_type
    
    case 'shuffle'

        [V,D] = eigenshuffle(A);
                
    case 'standard'
        
        D = NaN(n_dim, n_time);
        
        for tt = 1: n_time
            
            [~,dd] = eig(A(:,:,tt));
            dd = abs(diag(dd));
            [~,idx] = sort(dd,'descend');
            D(:,tt) = dd(idx);
    
        end     
end

rot_freq = angle(D)./(2*pi*dt);

end