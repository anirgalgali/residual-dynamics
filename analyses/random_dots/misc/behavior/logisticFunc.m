function [yy] = logisticFunc( xx, k, x0 )

yy = 1./(1 + exp(-k*(xx - x0))) ;


end