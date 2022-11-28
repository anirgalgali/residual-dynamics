function [data_out] = removeFsMultiplier(data)
%{ Sometimes the data is stored as firing rates instead of spike counts, 
%  and we would like to isntead work with spike counts. Since 
%  firing rate = spike_counts/dt = spike_counts*F, this function scales the
%  dataset by a factor (=Fs) if the data is indeed represented as rates and
%  not spike counts. 
% Input
% - data(struct) - containing  the neural data and associated variables for
%   a session, typically represented as binned temporall-varying firing rates
% Output
% - data_out(struct) - data scaled down by a factor of Fs.
% Author : Aniruddh Galgali (Jan 2018)
%}

% check if binned spike counts are in the right range (i.e each spike
% count is not multiplied by a multiplier (usually sampling freuency
% (Fs) is the multiplier)
    
data_out = data;
spCount_un = unique(data.response(:,:));
Fs = round(1./(data.time(2) - data.time(1)));
spCount_un_nz = spCount_un(spCount_un~=0);
if(mod(spCount_un_nz(1),Fs) == 0)

    data_out.response = (data.response)./Fs;
    
else
    return
end

end