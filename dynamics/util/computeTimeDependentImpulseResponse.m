function [out] = computeTimeDependentImpulseResponse(A,time,dirs,tRefs,varargin)

%{ This function computes a time-dependent impulse response of the dynamics
%  matrices specified in A
% Inputs
% - A (matrix/tensor) - either of size n_dyn_dim x n_dyn_dim
% (time-invariant) or n_dyn_dim x n_dyn_dim x n_times (time-varying)
% - time (array) - indicating the times associated with different A's
% - dirs (matrix) - of size ndyn_dim x n_dirs. Each column specifies a
% direction on state-space along which you provide the initial pulse for
% computing the impulse response.
% - tRefs (struct) - contains information about the times starting at which 
%   you want to compute the impulse response. 
% - varargin
% Outputs
% out (struct) - impulse response for different settings of tRefs/dirs.
% Author : Aniruddh Galgali (Written : Apr 2019)
%}
 
assert(size(A,3) == length(time),'# of dynamics matrices do not match number of time points');

numDirs = size(dirs,2);

time_labels = tRefs.time_labels(1) .* ones(length(time),1);
assignopts(who, varargin);

unique_time_labels = unique(time_labels);
tRef_count = 1;

for ii = 1: length(unique_time_labels)

    

    idx_times = time_labels == unique_time_labels(ii);

    A_seg = A(:,:,idx_times);

    time_seg = time(idx_times);

    

    tRefs_seg = tRefs.time(tRefs.time_labels == unique_time_labels(ii));
    tRefs_seg_labels = tRefs.time_labels(tRefs.time_labels == unique_time_labels(ii));
    n_times_ref = length(tRefs_seg);

    

    for jj = 1:n_times_ref

        

        [~,idx] = min(abs(time_seg - tRefs_seg(jj)));

        timeMask = false(1,length(time_seg));

        timeMask(idx:length(time_seg)) = true;

        dynNorm = zeros(numDirs,sum(timeMask));
        dynAngle = zeros(numDirs,sum(timeMask));
        propDynamics = zeros(size(A_seg,1),numDirs,sum(timeMask));        

        A_compose = eye(size(A_seg,1));

        tCount = 1;

        for tt = idx: length(time_seg)


            A_compose = A_seg(:,:,tt)*A_compose;

            

            for iDir = 1:numDirs

                

                dynamics_fwd_pass = A_compose * dirs(:,iDir);

        

                dynNorm(iDir,tCount) = norm(dynamics_fwd_pass);

                dynAngle(iDir,tCount) = normc(dynamics_fwd_pass)' * normc(dirs(:,iDir));

        

                propDynamics(:,iDir,tCount) = dynamics_fwd_pass;

        

            end

            

            tCount = tCount + 1;

    

        end

        

        out(tRef_count).norm = dynNorm;

        out(tRef_count).angle = dynAngle;

        out(tRef_count).prop_v = propDynamics;
        
        out(tRef_count).t_start = tRefs_seg(jj);

        out(tRef_count).t_start_label = tRefs_seg_labels(jj);
        
        time_axes = time_seg(timeMask);
        
        out(tRef_count).time = time_axes;

   

        tRef_count = tRef_count + 1;

    end


end

 

end