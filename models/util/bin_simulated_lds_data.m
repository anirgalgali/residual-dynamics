function [residuals_binned, data_binned] = bin_simulated_lds_data(X, time, bin_sizes, varargin)

if(nargin == 4)
    
    burn_in_size = varargin{1};

elseif(nargin == 3)

    burn_in_size = 0;

end
dt = time(2) - time(1);
data_binned = cell(numel(bin_sizes),1);
residuals_binned = cell(numel(bin_sizes),1);

for iB = 1:length(bin_sizes)
    
    data_S = [];
    data_S.response = X(:,burn_in_size + 1:end,:);
    data_S.time = time(burn_in_size + 1 : end);
    data_S = binSpikeCounts(data_S,bin_sizes(iB),1./dt,false,false);
    data_binned{iB} = data_S;
    residuals_binned{iB} = data_binned{iB}.response - repmat(mean(data_binned{iB}.response,3),[1 1 size(data_binned{iB}.response,3)]);
    
end
end