function data_smooth = tdrTemporalSmoothing(data,tmppars)
% Temporal smoothing of trial-by-trial or condition-averaged responses
%
% Inputs:
%  data: population response (simultaneous, see 'help tdr')
%  tmppars.filter: filter type ('gauss' or 'box')
%  tmppars.width: filter widht (s)
%
% Outputs:
%  data_smooth: smoothed data, same format as data.
%
% data_smooth = tdrTemporalSmoothing(data,tmppars)

% Initialize
data_smooth = data;

% Smooth
if isfield(tmppars,'width') && ~isempty(tmppars.width) && tmppars.width~=0
    % Sample duration
    dt = data.time(2) - data.time(1);
    
    % Make filter
    switch tmppars.filter
        case 'gauss'
            
            % Filter parameters
            fw = tmppars.width;
            
            % Time axis
            nt = round(fw*3.75/dt);  % this means we go out 3.75 times the standard deviation
            tf = (-nt:nt)*dt;
            
            % Make filter
            rf = exp(-tf.^2/(2*fw^2));
            rf = rf/(sum(rf)*dt);
            
        case 'box'
            
            % Filter parameters
            fw = tmppars.width;
            
            % Box width
            nt = round(fw/2/dt);
            
            % Make filter
            rf = ones(1,2*nt+1);
            rf = rf/(sum(rf)*dt);
    end
    
    % Temporal smoothing
    if length(rf)>1
        % Dimensions
        [nun npt ntr] = size(data.response);
        
        % Loop over units
        for iun = 1:nun
            % Loop over trials
            for itr = 1:ntr
                % Raw response
                rraw = squeeze(data.response(iun,:,itr));
                
                % Pad extremes
                rpad = [ones(1,length(rf))*rraw(1) rraw ones(1,length(rf))*rraw(end)];
                
                % Filter
                rlng = filter(rf,1,rpad)*dt;
                
                % Shift
                rfil = rlng(length(rf)+floor(length(rf)/2)+1:length(rf)+floor(length(rf)/2)+length(rraw));
                
                % Keep
                data_smooth.response(iun,:,itr) = rfil;
            end
        end
    end
end


