function vTime = tdrPassageOfTime(data,potpars)

% Dimensions
[ndm,npt,ntr] = size(data.response);

% Trials to use
if isfield(potpars,'trial_pot') && ~isempty(potpars.trial_pot)
    itrial = potpars.trial_pot;
else
    itrial = 1:ntr;
end

% Average response across all relevant trials
ravg = nanmean(data.response(:,:,itrial),3);

% Time window for urgency
[~,it1] = min(abs(potpars.time_win(1)-data.time));
[~,it2] = min(abs(potpars.time_win(2)-data.time));

% Time window for arc
itt = potpars.time_arc;

% Urgency direction
dir_urgency = ravg(:,it2)-ravg(:,it1);

% Arc direction
dir_arc = mean(ravg(:,itt) - repmat(((ravg(:,it2)-ravg(:,it1))/2 + ravg(:,it1)),[1 sum(itt)]),2);

% Define normalize vectors
response = zeros(ndm,1,2);
response(:,1,1) = normify(dir_urgency);
response(:,1,2) = normify(dir_arc);

% Keep what you need
vTime.name = {'urgency';'arc'};
vTime.response = response;
vTime.time = [];
vTime.dimension = data.dimension;

