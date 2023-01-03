function [Q,varargout] = determineLatentNoiseCovariance(C_X,A)

assert(size(C_X,1) == size(A,1),'wrong dimensions');
assert(size(C_X,3) == size(A,3) + 1,'time lengths do not match');
Q = NaN(size(A));
Q_org = NaN(size(A));
traceQ = NaN(size(A,3),1);
dd_org = NaN(size(A,1),size(A,3));
dd_near = NaN(size(A,1),size(A,3));

for tt = 1:size(A,3)
    
    Q_est = C_X(:,:,tt+1) - A(:,:,tt)*C_X(:,:,tt)*A(:,:,tt)' ;
    [vv,dd] = eig(Q_est);
    dd_org(:,tt) = sort(diag(dd),'descend');
    
    Q_org(:,:,tt) = Q_est;
    Q(:,:,tt) = mynearestSPD(Q_est);
    
    [vv,dd] = eig(Q(:,:,tt));
    dd_near(:,tt) = sort(diag(dd),'descend');
    
    traceQ(tt) = trace(Q(:,:,tt));
end
   
    
varargout{1} = traceQ;
varargout{2} = dd_org; 
varargout{3} = dd_near;
varargout{4} = Q_org;
end