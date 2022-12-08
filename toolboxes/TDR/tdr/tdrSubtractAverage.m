function data_sub = tdrSubtractAverage(data_raw)
% tdrSubtractAverage subtract average
%
% Input: 
%    data_raw: raw data
%
% Output:
%    data_sub: average-subtracted response
%
% data_sub = tdrSubtractAverage(data_raw)

% Initialize
data_sub = data_raw;

% Raw responses
Rraw = data_raw.response;

% Dimensions
[ndm,npt,ncd] = size(Rraw);

% Average response
Ravg = repmat(nanmean(Rraw,3),[1 1 ncd]);

% Average-subtracted response
data_sub.response = Rraw - Ravg;