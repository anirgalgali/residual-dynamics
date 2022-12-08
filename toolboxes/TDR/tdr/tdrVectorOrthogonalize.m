function [coef_ort,lowUN_lowTA] = tdrVectorOrthogonalize(coef_fix,ortpars)
% tdrVectorOrthogonalize orthogonalize regression vectors
%
% Inputs:
%  coef_fix: regression vectors
%     .name: vector names {nvc 1}
%     .response: vector coefficients [ndm npt nvc]
%     .dimension: dimension labels {ndm 1}
%  ortpars.name: vectors to orthogonalize {nrg 1}
% 
% Outputs:
%  coef_ort: orthogonalized regression vectors
%     .name: vector names {nrg 1}
%     .response: vector coefficients [ndm npt nrg]
%     .dimension: dimension labels {ndm 1}
%  lowUN_lowTA: projection matrix from task-related subspace (subspace
%     basis) into original state space (original basis)
%  
% [coef_ort,lowUN_lowTA] = tdrVectorOrthogonalize(coef_fix,ortpars)

% Dimensions
nrg = length(ortpars.name);
[nun,npt,~] = size(coef_fix.response);

% Initialize
coef_ort = coef_fix;
coef_ort.name = ortpars.name;
coef_ort.response = zeros(nun,npt,nrg);

% Loop over time points
for ipt = 1:npt
    raw = zeros(nun,nrg);
    for irg = 1:nrg
        % Find vector
        jmatch = strcmp(coef_fix.name,ortpars.name{irg});
        
        % Keep vector
        raw(:,irg) = coef_fix.response(:,ipt,jmatch);
    end
    
    % Orthogonalize
    [qq,rr] = myqr(raw);
    ort = qq(:,1:nrg);
    
    % Keep what you need
    coef_ort.response(:,ipt,:) = ort;
end

% Keep the projection matrices
lowUN_lowTA = coef_ort.response;
