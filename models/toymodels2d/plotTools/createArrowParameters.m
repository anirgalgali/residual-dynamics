function [arrow_shape_scaled] = createArrowParameters(meta_params,l)
%{ This function creates a set of parameters that determine the graphical 
% charasterics of arrows that are ultimately used to plot the flow-fields.
% Play around with the following 8 hyper-parameters to create 
% aesthetically pleasing arrows/arrowheads. The hyper-parameters control 
% the different aspects of the arrow, such as tail length, 
% head to tail ratio, width of arrow head, etc..
%
%
% meta_params(1) = k_w
% meta_params(2) = k_j
% meta_params(3) = l1;
% meta_params(4) = l2;
% meta_params(5) = a1;
% meta_params(6) = a2
% meta_params(7) = s1;
% meta_params(8) = s2
%
% Author: Aniruddh Galgali (Oct 2018)

if( l <= meta_params(3))
    

    a_bar = meta_params(5);
 

elseif(l >= meta_params(4))
    

    a_bar = meta_params(6);
 

else
    

    slope_a = (meta_params(6) - meta_params(5))./(meta_params(4) - meta_params(3));
    a_bar = slope_a*(l-meta_params(3)) + meta_params(5);
 

end
 

if( l <= meta_params(3))
    

    s_bar = meta_params(7);
    

elseif(l >= meta_params(4))
    

    s_bar = meta_params(8);
    

else
    

    slope_s = (meta_params(8) - meta_params(7))./(meta_params(4) - meta_params(3)); 
    s_bar = slope_s*(l-meta_params(3)) + meta_params(7);
 

end
 

a = a_bar./l;
s = s_bar./l;
headW = meta_params(1).*a;
headL = a;
headI = meta_params(2)*a;
lineW = s;
 

arrow_shape_scaled = [headW headL headI lineW];

end