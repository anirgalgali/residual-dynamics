function [A,varargout] = generate_strongly_nonnormal_matrix(n_dim, tau_min, tau_max, dt, do_continuous)

ev_min = exp(-dt/tau_min);
ev_max = exp(-dt/tau_max);
s1 = 0;
angle_with_largest = 0;
count = 0;

while(s1 <= 1  |  abs(angle_with_largest - 90) > 15)

    V = normc(randn(n_dim));
    D = ev_min + (ev_max - ev_min) * rand(n_dim,1);
    Dsort = sort(D,'descend');
    A = V * diag(D) * inv(V);
    [uu,ss,~] = svd(A);
    s1 = ss(1,1);
    s2 = ss(2,2);

    [~,ind] = sort(D,'descend');
    angle_with_v = rad2deg(acos(uu(:,1)'* V(:,ind)));
    angle_with_largest = angle_with_v(1);
    count = count+ 1;
    
    if(s2 > 1 || s1 - 1 < 0.1 || Dsort(1) - Dsort(2) < 0.05)
        s1 = 0;
        angle_with_largest = 0;
        continue
    end
    
    if(count > 10000)
        fprintf('search failed\n')
        A = [];
        break
    end
end

if(do_continuous)
   if(~isempty(A))
       
       A = logm(A)./dt;
       
   end
end
a_stats.s1 = s1;
a_stats.s2 = s2;
a_stats.angle_llsv_ev = angle_with_v;
a_stats.D = D;

varargout{1} = a_stats;
varargout{2} = count;
end