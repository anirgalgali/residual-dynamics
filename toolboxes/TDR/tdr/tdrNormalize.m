function data_nrml = tdrNormalize(data,nrmlpars)
% tdrNormalize Normalize population responses
%
% Computes R* = (R - m)/(s + c)
% 
%  R* are the normalized responses
%  R are the raw responses
%  m is the mean across trials/conditions and time
%  s is the standart deviation across trials/conditions and time
%  c is a constant, to prevent R* from exploding for s->0
%
% Inputs:
%  data: population response (sequential or simultaneous, see 'help tdr')
%  nrmlpars.ravg: the mean across trials/conditions and time [nun 1]
%  nrmlpars.rstd: the std across trials/conditions and time [nun 1]
%  nrmlpars.cnst: the constant [nun 1] or [1 1] (def = 0)
%
% Ouputs:
%  data_nrml: normalized responses
%
% data_nrml = tdrNormalize(data,nrmlpars)


if nargin < 2
    nrmlpars = [];
end

% Check if data is sequentially or simultaneously recorded
if isfield(data,'unit') && ~isempty(data.unit)
    
    
    %--- Sequential recordings ---    
    % Number of units
    nun = length(data.unit);
    
    % Initialize
    data_nrml = data;
    
    % Parameters
    if isempty(nrmlpars) || ~isfield(nrmlpars,'ravg') || isempty(nrmlpars.ravg)
        nrmlpars.ravg = zeros(nun,1);
    end
    if isempty(nrmlpars) || ~isfield(nrmlpars,'rstd') || isempty(nrmlpars.rstd)
        nrmlpars.rstd = ones(nun,1);
    end
    if isempty(nrmlpars) || ~isfield(nrmlpars,'cnst') || isempty(nrmlpars.cnst)
        nrmlpars.cnst = zeros(nun,1);
    elseif length(nrmlpars.cnst)==1
        nrmlpars.cnst = ones(nun,1)*nrmlpars.cnst;
    end
        
    % Loop over units
    for iun = 1:nun
        
        % Z-score
        data_nrml.unit(iun).response = ...
            (data.unit(iun).response - nrmlpars.ravg(iun)) / ...
            (nrmlpars.rstd(iun) + nrmlpars.cnst(iun));
        
    end
    
else
    
    
    %--- Simultaneous recordings ---
    % Initialize
    data_nrml = data;

    % Dimensions
    [nun npt ncd] = size(data.response);
    
    % Parameters
    if isempty(nrmlpars) || ~isfield(nrmlpars,'ravg') || isempty(nrmlpars.ravg)
        nrmlpars.ravg = zeros(nun,1);
    end
    if isempty(nrmlpars) || ~isfield(nrmlpars,'rstd') || isempty(nrmlpars.rstd)
        nrmlpars.rstd = ones(nun,1);
    end
    if isempty(nrmlpars) || ~isfield(nrmlpars,'cnst') || isempty(nrmlpars.cnst)
        nrmlpars.cnst = zeros(nun,1);
    elseif length(nrmlpars.cnst)==1
        nrmlpars.cnst = ones(nun,1)*nrmlpars.cnst;
    end
    
    % Mean and STD
    Ravg = repmat(nrmlpars.ravg,[1 npt ncd]);
    Rstd = repmat(nrmlpars.rstd,[1 npt ncd]);
    Cnst = repmat(nrmlpars.cnst,[1 npt ncd]);
    
    % Z-score
    data_nrml.response = (data.response - Ravg) ./ (Rstd + Cnst);
    
end

