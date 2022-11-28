function [data_out] = applyUnitRemoval(data,unitMask)
%{ Removes units/dimensions from dataset based on units flagged in unitMask
% Input
% - data(struct) - containing responses and associated task/trial variables
% - unitMask (array, logical) - of size nUnits x 1, with 1s indicating
%   units to be removed
% Output
% - data_out (struct) - filtered data struct with removed units
%
% Author : Aniruddh Galgali (May 2017)
%}
data_out = data;
data_out.response = [];
data_out.dimension = [];
data_out.response = data.response(~unitMask,:,:);
data_out.dimension = data.dimension(~unitMask);
if(isfield(data_out,'dspInfo'))
    data_out.dspInfo = data.dspInfo(~unitMask);
end