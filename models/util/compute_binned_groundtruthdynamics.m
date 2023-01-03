function [binned_A, varargout] = compute_binned_groundtruthdynamics(A, step_size, bin_size)

nBins = round(bin_size/(step_size * 1000));

if(ndims(A) == 3)
    
    n_times = size(A,3);
    time = 0 : step_size : (n_times - 1)*step_size;
    T_bins    = floor(size(A,3)/ nBins);
    binned_A = NaN(size(A,1),size(A,1),T_bins);
    binned_time = NaN(1,T_bins);

    % This bins the ground truth matrices that have been defined at a finer
    % temporal resolution
    for tt = 1:T_bins

        iStart = nBins * (tt-1) + 1;
        iEnd   = nBins * tt;
        avg_A =  mean(A(:, :, iStart:iEnd), 3);
        binned_A(:,:,tt) = discretize_dynamics_matrix(avg_A, step_size * bin_size);
        binned_time(tt) = median(time(iStart:iEnd));
    end
    
    varargout{1} = binned_time;
    
elseif(ndims(A) == 2)
    
    binned_A = discretize_dynamics_matrix(A, step_size * bin_size);
    
    
end


end