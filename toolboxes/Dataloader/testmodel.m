function rr = testmodel(pars)
% a fake model 

switch (nargin)
	case(9)
    case(0)
        rr = {'\sigma_c','Rmax','\sigma_s/\sigma_c','g','c50','\sigma_i','kmask','\sigma_f'};
        return;
    otherwise
        error('Need to provide 0 or 9 args')
end