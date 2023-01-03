function [dd_out, varargout] =  rebin_eigenvalue_dynamics(A, dt_in, dt_out,sort_type)

%units of dt_out and dt_in are in seconds
% A should be a discrete time dynamics matrix

dd_out = NaN(size(A,1),size(A,3));
% [~,Dseq] = eigenshuffle(A);
switch sort_type
    
    case 'standard'
        
        for tt =1:size(A,3)
            
            [~,dd_in] = eig(A(:,:,tt));
            [~, idx] = sort(abs(diag(dd_in)),'descend');
            dd_in = diag(dd_in);
            dd_in = dd_in(idx);
            dd_out(:,tt) = dd_in.^(dt_out/dt_in);
            
        end
        
    case 'shuffle'
        
        [Vseq,Vseq_r,Dseq] = my_eigenshuffle(A,'abs');
        [Useq,Sseq,Vseq] = my_singularshuffle(A,  'right');
        dd_out = Dseq.^(dt_out/dt_in);
        ss_out = Sseq.^(dt_out/dt_in);
        varargout{1} = ss_out;
end

% 
end