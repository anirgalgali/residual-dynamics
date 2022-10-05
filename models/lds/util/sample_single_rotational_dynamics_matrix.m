function [A] = sample_single_rotational_dynamics_matrix(tau, rot_f, vy_re_max, vy_re_min, vy_im_max, vy_im_min, dt, varargin)

if(nargin == 8)
    rseed = varargin{1};
    rng(rseed);
end

r = 2*pi*rot_f*dt;
g = exp(-dt/tau);
a = g*cos(r);
b = g*sin(r);
ev = [a + 1i*b; a - 1i*b];
c = vy_re_min + (vy_re_max - vy_re_min)*rand(1);
d = vy_im_min + (vy_im_max - vy_im_min)*rand(1);
eig_vec = [1 1; c + 1i*d c - 1i*d];
A = eig_vec * diag(ev) * inv(eig_vec);

end


