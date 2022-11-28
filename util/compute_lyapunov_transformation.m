function [ss] = compute_lyapunov_transformation(A)
%dt is in ms
[~,p,T] = size(A);
% converting to continuous time representation
phiA = eye(p);
ss = NaN(p,T+1);
c_count = 1;
for tt = 1:T
    
    phiA = A(:,:,tt)*phiA;
    [qq,rr] = myqr(phiA);
    D = diag(sign(diag(rr)));
    qq = qq*D;
    rr = D*rr;
    
    phiA = qq; 
    ss(:,tt) = log(diag(rr));
    c_count = c_count  + 1;
end


end
