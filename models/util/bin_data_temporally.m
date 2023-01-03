function [data_binned] = bin_data_temporally(data, bin_sizes, varargin)

if(nargin == 3)
    
    burn_in_size = varargin{1};

elseif(nargin == 2)

    burn_in_size = 0;

end

data_binned = repmat(data,[1 length(bin_sizes)]);
dt = data.time(2) - data.time(1);

for iB = 1:length(bin_sizes)
    
    data_S = [];
    data_S.response = data.response(:,burn_in_size + 1:end,:);
    data_S.time = data.time(burn_in_size + 1 : end);
    data_S.time_iev = data.time_iev(burn_in_size + 1 : end);
    data_S.time_rel = data.time_rel(burn_in_size + 1 : end);
    
    data_S = binSpikeCounts(data_S,bin_sizes(iB),1./dt,false,false);
    
    data_binned(iB).response = data_S.response;
    data_binned(iB).time = data_S.time;
    data_binned(iB).time_iev = data_S.time_iev;
    data_binned(iB).time_rel = data_S.time_rel;
   
end


end