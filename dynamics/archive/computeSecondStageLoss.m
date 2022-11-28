function [loss] = computeSecondStageLoss(Y,X,A)
    
[p,T,K] = size(Y);

D = spdiags(ones(p.^2*(T-1),1)*[-1 1],[0 p.^2],p.^2*(T-1),p.^2*(T));
Y_v = NaN(p*K*T,1);
for tt = 1:size(X,2)
    
    B{tt} = sparse(kron(squeeze(X(:,tt,:))', eye(p)));
    Y_v((tt-1)*p*K + 1 : tt*p*K) = vec(Y(:,tt,:));
    
end


Ztilde = blkdiag(B{:});
loss = norm(Y_v(:) - Ztilde*vec(A)).^2 ;

end