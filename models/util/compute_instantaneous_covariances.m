function [covX, totalCovX] =  compute_instantaneous_covariances(X)

X = squeeze(num2cell(permute(X,[3 1 2]),[1 2]));
covX = cellfun(@(x) cov(x),X,'uni',false);
totalCovX = cell2mat(cellfun(@(x) trace(x),covX,'uni',false));
covX = cat(3, covX{:});      
end