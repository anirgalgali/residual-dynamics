function [l0,l1,r,p_var,varargout] = compute_prop_ratio(X)

[~,n_t,n_tr] = size(X);
s0_all = zeros(n_t,1);
s1_all = zeros(n_t,1);
prop_var = zeros(n_t,1);

if(size(X,1)>1)
    for ii = 1:n_t
        
        C0 = (1./n_tr).*(squeeze(X(:,ii,:))*squeeze(X(:,ii,:))');
        if(ii < n_t)
            C1 = (1./n_tr).*(squeeze(X(:,ii+1,:))*squeeze(X(:,ii,:))');
        end
        
        [u0,s0,v0] = svd(C0);
        [u1,s1,v1] = svd(C1);
        
        s0_all(ii) = sum(diag(s0).^2);
        s1_all(ii) = sum(diag(s1).^2);
        prop_var(ii) = s1(1,1);
    end
    
    l0 = mean(s0_all);
    l1 = mean(s1_all);
    p_var = mean(prop_var);
    
    varargout{1} = s0_all;
    varargout{2} = s1_all;
    varargout{3} = prop_var;
    varargout{4} = s1_all./s0_all;


elseif(size(X,1) == 1)
   
    t_cov = cov(squeeze(X)');
    C0 = diag(t_cov);
    C1 = diag(t_cov,1);
    
    l0 = mean(C0(2:end));
    l1 = mean(C1);
    p_var = l1;

end


r = l1./l0;



end