function [Summary] = computejPCA(X, time_rel, trial_labels, time_labels, analyze_time_mask, n_jpc_planes, varargin)
%{ This function computes the jPC planes given some neural data. The jPC planes
%  capture rotational dynamics in neural trajectories.
%  Input
%  - X (matrix) - of single trial data of size n_dim x n_times x n_trials.
%   Typically, we stitch responses from different epochs (if any) along the
%   time dimension, therefore X contains responses from both epochs. However,
%   jPCA is computed separately for each epoch.
%  -time_rel (array) - of size 1 x n_times,  indicating the 'relative' times
%   in each epoch. 
%  -trial_labels (array) - of size 1 x n_trials,  indicating the task condition
%   associated with each trial
%  -time_labels (array) - of size 1 x n_times,  indicating the epoch of
%   each time-point.
%  -analyze_time_mask (cell array) of size n_align x 1 , where n_align is
%   the total number of alignments/task epochs (i.e number of unique
%   elements of time_labels). Each cell array consists of a time -window
%   specified as [tStart tEnd], which indicates the range of the time window
%   to consider while computing jPCs
%  - n_jpc_panes - specifies the number of JPC planes to extract
%  - varargin
%     - filt_pars (struct) - structure containing additional information
%     about the filter used to smooth responses prior to JPC computation.
%
%  Output
%  - SUmmary(struct) - contains all the information about the JPCs and the
%    respective projections
%  
%  This code is an adapation of the original jPCA toolbox feveloped by Mark 
%  Churchland/ John Cunnningham 
%  (see https://churchland.zuckermaninstitute.columbia.edu/content/code)
%
%  Modified: ANiruddh Galgali (Jan 2021)
%}
nPCs = 2*n_jpc_planes; 

unique_trial_labels= unique(trial_labels);
nconds = length(unique_trial_labels);

unique_time_labels = unique(time_labels);
nalign = length(unique_time_labels);

if(nargin == 7)
    
    filt_pars = varargin{1};
    assert(isfield(filt_pars,'filt_width'),'need to specify width of smoothing filter');
    data_S = [];
    
    for ialign = 1:nalign
        
        time_ev = time_rel(time_labels == unique_time_labels(ialign)); 
        data_S.response = X(:,time_labels == unique_time_labels(ialign),:);
        data_S.time = time_ev;
        data_S = smoothResponses(data_S,filt_pars.filt_width,time_ev(2) - time_ev(1));
        X(:,time_labels == unique_time_labels(ialign),:) = data_S.response ;   
        
    end
    
end

assert(length(analyze_time_mask) == nalign);
Summary = repmat(struct(),[nalign 1]);

for ialign = 1:nalign
    analyzeIndices = analyze_time_mask{ialign};
    if size(analyze_time_mask{ialign},1) == 1
        analyzeIndices = analyzeIndices';  % orientation matters for the repmat below
    end
    
    analyzeMask = repmat(analyzeIndices,2,1); 
    bunchOtruth = true(sum(analyzeIndices)-1,1);
    maskT1 = repmat( [bunchOtruth;false],2,1);  % skip the last time for each condition
    maskT2 = repmat( [false;bunchOtruth],2,1);  % skip the first time for each condition
    
    time_ev = time_rel(time_labels == unique_time_labels(ialign));
    
    if sum(analyzeIndices) < 5
        disp('warning, analyzing few or no times');
        
    end
    
    X_ev = X(:,time_labels == unique_time_labels(ialign),:);
    X_ev_avg = NaN(size(X_ev,1), size(X_ev,2),nconds);
    for icond = 1:nconds
        
       X_ev_avg(:,:,unique_trial_labels(icond)) = mean(X_ev(:,:,trial_labels == unique_trial_labels(icond)),3);
        
    end
    
    bigA =  X_ev_avg(:,:)';
    smallA = bigA(analyzeMask,:);
    [PCvectors,rawScores] = pca(smallA);
    
    PCvectors = PCvectors(:, 1 : nPCs);  % Want to always find atleast  jPC planes 
    Ared = rawScores(:,1:nPCs);
    
    numAnalyzedTimes = size(Ared,1)/2;
    dState = Ared(maskT2,:) - Ared(maskT1,:);  % the masks just give us earlier and later times within each condition
    dState = dState./(time_ev(2) - time_ev(1)); % scaling by dt
    preState = Ared(maskT1,:);  % just for convenience, keep the earlier time in its own variable
    
    M = (dState'/preState');  % M takes the state and provides a fit to dState
    
    Mskew = skewSymRegress(dState,preState)';  % this is the best Mskew for the same equation
    
    [V,D] = eig(Mskew); % V are the eigenvectors, D contains the eigenvalues
    evals = diag(D); % eigenvalues
    [~,sortIndices] = sort(abs(evals),1,'descend');
    evals = evals(sortIndices);  % reorder the eigenvalues
    evals = imag(evals);  % get rid of any tiny real part
    V = V(:,sortIndices);  % reorder the eigenvectors (base on eigenvalue size)
    
    jPCs = zeros(size(V));
    
    for pair = 1:ceil(size(Ared,2)/2)
        
        vi1 = 1+2*(pair-1);
        vi2 = 2*pair;
        
        VconjPair = V(:,[vi1,vi2]);  % a conjugate pair of eigenvectors
        evConjPair = evals([vi1,vi2]); % and their eigenvalues
        VconjPair = getRealVs(VconjPair, evConjPair, Ared, numAnalyzedTimes);
        
        jPCs(:,[vi1,vi2]) = VconjPair;
    end
    
    bigAred = bsxfun(@minus, bigA, mean(bigA)) * PCvectors(:,1:nPCs);

    origVar = trace(cov(bigA));
    varCaptEachPC = diag(PCvectors(:,1:nPCs)' * cov(bigA) * PCvectors(:,1:nPCs))./ origVar ;
    varCaptEachJPC = diag((PCvectors(:,1:nPCs)*jPCs)' * cov(bigA) * (PCvectors(:,1:nPCs)*jPCs))./ origVar ;
    varCaptEachPlane = reshape(varCaptEachJPC, 2, nPCs/2);
    varCaptEachPlane = sum(varCaptEachPlane);
    
    
    numAnalyzedTimes = size(Ared,1)/nconds;
    
    projAllTimes = bigAred * jPCs;
    tradPCA_AllTimes = bsxfun(@minus, bigA, mean(smallA)) * PCvectors;  % mean center in exactly the same way as for the shorter time period.

    index1 = 1;
    index2 = 1;
    Projection = repmat(struct(),[nconds 1]);
    
    for c = 1:nconds
        
        index1b = index1 + numAnalyzedTimes -1;  % we will go from index1 to this point
        index2b = index2 + length(time_ev) -1;  % we will go from index2 to this point
        Projection(c).times = time_ev(analyzeIndices);
        Projection(c).projAllTimes = projAllTimes(index2:index2b,:);
        Projection(c).allTimes = time_ev;
        Projection(c).tradPCAproj = Ared(index1:index1b,:);
        Projection(c).tradPCAprojAllTimes = tradPCA_AllTimes(index2:index2b,:);
        index1 = index1+numAnalyzedTimes;
        index2 = index2+length(time_ev);
    
    end


   
    Summary(ialign).jPCs = jPCs;
    Summary(ialign).evals = evals;
    Summary(ialign).PCs = PCvectors;
    Summary(ialign).jPCs_fullD = PCvectors * jPCs;
    Summary(ialign).varCaptEachJPC = varCaptEachJPC;
    Summary(ialign).varCaptEachPC = varCaptEachPC;
    Summary(ialign).varCaptEachPlane = varCaptEachPlane;
    Summary(ialign).Mbest = M;
    Summary(ialign).Mskew = Mskew;
    Summary(ialign).Projection = Projection;   
end

end

function Vr = getRealVs(V,evals,Ared,numAnalyzedTimes)

    % get real vectors made from the eigenvectors
    
    % by paying attention to this order, things will always rotate CCW
    if abs(evals(1))>0  % if the eigenvalue with negative imaginary component comes first
        Vr = [V(:,1) + V(:,2), (V(:,1) - V(:,2))*1i]; 
    else
        Vr = [V(:,2) + V(:,1), (V(:,2) - V(:,1))*1i];
    end
    Vr = Vr / sqrt(2);

    % now get axes aligned so that plan is spread mostly along the horizontal axis
    testProj = (Vr'*Ared(1:numAnalyzedTimes:end,:)')'; % just picks out the plan times
    rotV = pca(testProj,'Economy',false);
    crossProd = cross([rotV(:,1);0], [rotV(:,2);0]);
    if crossProd(3) < 0, rotV(:,2) = -rotV(:,2); end   % make sure the second vector is 90 degrees clockwise from the first
    Vr = Vr*rotV; 

    % flip both axes if necessary so that the maximum move excursion is in the positive direction
    testProj = (Vr'*Ared')';  % all the times
    if max(abs(testProj(:,2))) > max(testProj(:,2))  % 2nd column is the putative 'muscle potent' direction.
        Vr = -Vr;
    end
end

function phaseData = getPhase(Proj, whichPair)
    numConds = length(Proj);
    d1 = 1 + 2*(whichPair-1);
    d2 = d1+1;
    
    for c=1:numConds
        data = Proj(c).proj(:,[d1,d2]);
        phase = atan2(data(:,2), data(:,1));  % Y comes first for atan2
        
        deltaData = diff(data);
        phaseOfDelta = atan2(deltaData(:,2), deltaData(:,1));  % Y comes first for atan2
        phaseOfDelta = [phaseOfDelta(1); phaseOfDelta];  %#ok<AGROW> % so same length as phase
        radius = sum(data.^2,2).^0.5;
        
        % collect and format
        % make things run horizontally so they can be easily concatenated.
        phaseData(c).phase = phase'; %#ok<AGROW>
        phaseData(c).phaseOfDelta = phaseOfDelta'; %#ok<AGROW>
        phaseData(c).radius = radius'; %#ok<AGROW>
        
        % angle between state vector and Dstate vector
        % between -pi and pi
        phaseData(c).phaseDiff = minusPi2Pi(phaseData(c).phaseOfDelta - phaseData(c).phase); %#ok<AGROW>
    end
    
end

%% plotting the phase difference between dx(t)/dt and x(t) where x is the 2D state
function circStatsSummary = plotPhaseDiff(phaseData, params, jPCplane)
    % compute the circular mean of the data, weighted by the r's
    circMn = circ_mean([phaseData.phaseDiff]', [phaseData.radius]');
    resultantVect = circ_r([phaseData.phaseDiff]', [phaseData.radius]');
    
    
    bins = pi*(-1:0.1:1);
    cnts = histc([phaseData.phaseDiff], bins);  % not for plotting, but for passing back out
    

    % do this unless params contains a field 'suppressHistograms' that is true
    if ~isfield(params,'suppressHistograms') || ~params.suppressHistograms
        figure;
        hist([phaseData.phaseDiff], bins); hold on;
        plot(circMn, 20, 'ro', 'markerFa', 'r', 'markerSiz', 8);
        plot(pi/2*[-1 1], [0 0], 'ko', 'markerFa', 'r', 'markerSiz', 8);
        set(gca,'XLim',pi*[-1 1]);
        title(sprintf('jPCs plane %d', jPCplane));
    end
    
    %fprintf('(pi/2 is %1.2f) The circular mean (weighted) is %1.2f\n', pi/2, circMn);
    
    % compute the average dot product of each datum (the angle difference for one time and condition)
    % with pi/2.  Will be one for perfect rotations, and zero for random data or expansions /
    % contractions.
    avgDP = averageDotProduct([phaseData.phaseDiff]', pi/2);
    %fprintf('the average dot product with pi/2 is %1.4f  <<---------------\n', avgDP);
    
    circStatsSummary.circMn = circMn;
    circStatsSummary.resultantVect = resultantVect;
    circStatsSummary.avgDPwithPiOver2 = avgDP;  % note this basically cant be <0 and definitely cant be >1
    circStatsSummary.DISTRIBUTION.bins = bins;
    circStatsSummary.DISTRIBUTION.binCenters = (bins(1:end-1) + bins(2:end))/2;
    circStatsSummary.DISTRIBUTION.counts = cnts(1:end-1);
    circStatsSummary.RAW.rawData = [phaseData.phaseDiff]';
    circStatsSummary.RAW.rawRadii = [phaseData.radius]';
    
end

%% for computing the average dot product with a comparison angle
function avgDP = averageDotProduct(angles, compAngle, varargin)
    
    x = cos(angles-compAngle);
    
    if ~isempty(varargin)
        avgDP = mean(x.*varargin{1}) / mean(varargin{1});  % weighted sum
    else
        avgDP = mean(x);
    end
end