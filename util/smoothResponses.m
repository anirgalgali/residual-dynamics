function [data_out] = smoothResponses(data, filtWidth, stepSize)
%{ Function to smooth neural spike counts/firing rates
% Input 
% - data (struct) - contains the neural data and the associated task
%  variable for a single recording session.
% - filtWidth(scalar) - specifies the width of the smoothing filter
%  (in seconds)
% - stepSize (scalar) - time interval between adjacent time bins in data
% Outputs
% - data_out (struct) - smoothed data output
%
% Author: Aniruddh Galgali. Heavily adapted from the TDR toolbox by Valerio
% Mante.
%}
  fltHL = ceil((filtWidth/2) / stepSize);
 
  flt = ones(1,2*fltHL+1);
  
  yDim = size(data.response,1);
  T = size(data.response,2);
  K = size(data.response,3);
  nm = conv(flt, ones(1, T)); 
  
  data_out = data;

  % Normalize by sum of filter taps actually used
  
for kk = 1:K  
  for ii = 1:yDim
    ys = conv(flt, squeeze(data.response(ii,:,kk))) ./ nm;
    % Cut off edges so that result of convolution is same length 
    % as original data
    data_out.response(ii,:,kk) = ys(fltHL+1:end-fltHL);
  end
end


%   data_smooth = data;
% 
% % Smooth
% if isfield(tmppars,'width') && ~isempty(tmppars.width) && tmppars.width~=0
%     % Sample duration
%     dt = data.time(2) - data.time(1);
%     
%     % Make filter
%     switch tmppars.filter
%         case 'gauss'
%             
%             % Filter parameters
%             fw = tmppars.width;
%             
%             % Time axis
%             nt = round(fw*3.75/dt);  % this means we go out 3.75 times the standard deviation
%             tf = (-nt:nt)*dt;
%             
%             % Make filter
%             rf = exp(-tf.^2/(2*fw^2));
%             rf = rf/(sum(rf)*dt);
%             
%         case 'box'
%             
%             % Filter parameters
%             fw = tmppars.width;
%             
%             % Box width
%             nt = round(fw/2/dt);
%             
%             % Make filter
%             rf = ones(1,2*nt+1);
%             rf = rf/(sum(rf)*dt);
%     end
%     
%     % Temporal smoothing
%     if length(rf)>1
%         % Dimensions
%         [nun npt ntr] = size(data.response);
%         
%         % Loop over units
%         for iun = 1:nun
%             % Loop over trials
%             for itr = 1:ntr
%                 % Raw response
%                 rraw = squeeze(data.response(iun,:,itr));
%                 
%                 % Pad extremes
%                 rpad = [ones(1,length(rf))*rraw(1) rraw ones(1,length(rf))*rraw(end)];
%                 
%                 % Filter
%                 rlng = filter(rf,1,rpad)*dt;
%                 
%                 % Shift
%                 rfil = rlng(length(rf)+floor(length(rf)/2)+1:length(rf)+floor(length(rf)/2)+length(rraw));
%                 
%                 % Keep
%                 data_smooth.response(iun,:,itr) = rfil;
%             end
%         end
%     end
end
