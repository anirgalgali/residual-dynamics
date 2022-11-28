function [sVal] = computeTimeVaryingHankelSingularValues(H)

if(ndims(H) == 3)
    
    for tt = 1:size(H,3)
        
        [~,S,~] = svd(H(:,:,tt),'econ');
        sVal(:,tt) = diag(S);
    end
    
elseif(ndims(H) == 2)
    
    [~,S,~] = svd(H,'econ');
    sVal= diag(S);
    
end

end