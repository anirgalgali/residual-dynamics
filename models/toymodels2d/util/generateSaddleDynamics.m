function [model_out] =  generateSaddleDynamics(pars,type)
%{ This function returns the dynamics update equation corresponding to 
% a saddle point model. The model can either be characterized by an
% unstable fixed point at the origin (simple_saddle) or an unstable fixed
% point in an arbitrary location (complex_saddle). The function also
% analytically computes the fixed points of the model.
%}
switch type
    
    case 'simple_saddle'
        
        syms x y
        x_dot = [(pars.a)*x^3 + (pars.b)*x; pars.c*y];
        model_ode = matlabFunction(x_dot, 'Vars', {[x, y]});
        J = [pars.b + 3*pars.a*x.^2 0;0 pars.c];
        model_jac = matlabFunction(J,'Vars', {[x, y]});
        
        
        r = roots([pars.a 0 pars.b 0]);
        
        idx_roots_valid = find(imag(r) == 0); % keep only real roots
        
        for ii = 1: length(idx_roots_valid)
            fp_out(ii).fp_loc(1,1) = r(idx_roots_valid(ii));
            fp_out(ii).fp_loc(2,1) = 0;
            [fp_out(ii).eVec,fp_out(ii).eVal] = eig(model_jac(fp_out(ii).fp_loc'));
        end
        
        
        
        
    case 'complex_saddle'
        
        
        syms x y
        x_dot = [-y + pars.b*x + pars.a*x.^3;-x + pars.c*y.^3];
        model_ode = matlabFunction(x_dot, 'Vars', {[x, y]});
        
        J = [pars.b + 3*pars.a*x.^2 -1;-1 3*pars.c*y.^2];
        model_jac = matlabFunction(J,'Vars', {[x, y]});
        
        
        r = roots([(pars.c)^3*pars.a 0 0 0 0 0 pars.b*pars.c 0 -1 0]);
        idx_roots_valid = find(imag(r) == 0); % keep only real roots
        for ii = 1: length(idx_roots_valid)
            
            fp_out(ii).fp_loc(1,1) = pars.c*r(idx_roots_valid(ii)).^3;
            fp_out(ii).fp_loc(2,1) = r(idx_roots_valid(ii));
            [fp_out(ii).eVec,fp_out(ii).eVal] = eig(model_jac(fp_out(ii).fp_loc'));
            
        end
        
end


model_out.model_ode = model_ode;
model_out.model_jac = model_jac;
model_out.fixed_points = fp_out;


end

