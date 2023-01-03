function [analysis] = computeTimeDependentFlowField_v2(A,U,V,gridX,gridY,avg_trajectory,single_trials,time_iev,time_rel,time_abs)
%{ This function computes a time-dpendent flow-field based on the estimate of residual dynamics 
% (specified in A) at different locations in a 2D grid.
% Input
% - A (3d tensor) - n_duyn_dim x n_dyn_dim x n_times-  estimates of
%   residual dynamics at each time. 
% - U (matrix) - specifies the dynamics subspace
% - V (matrix) - of size n_dim x 2, where each column of V specifies a
%   direcion in state-space. The two directions together constitute a plane 
%   in which we want to compute the flow-field.
% - gridX,gridY (3d tensor) - specifies locations on a grid in the subspace
%   spanned by columns of V. The first dimension specifies the x and y
%   coordinate, and the second/third dimensions indexes the different
%   locations.
% - avg_trajectory (3d tensor) - specifying the average trajectory
% - single_trials (3d tensor) - specifying the single trials
% - time_iev (array) - label indicating the temporal alignment that each
%   time point belongs to 
% - time_rel (array) - times relative to each temporal alignment
% - time_abs (array) - absolute times.
%
% Output
% - analysis (struct) - containing the estimates of the flow-field.
%} 
assert(iscell(A),'dynamics needs to be passed as a cell array')
assert(iscell(V),'projection planes need to be passed as a cell array');

numChoices = length(A);
assert(numChoices == size(avg_trajectory,3),'number of average trajectories should match');

[grid_locs_X,grid_locs_Y] = meshgrid(gridX,gridY);
grid_points = cat(3,grid_locs_X,grid_locs_Y);
grid_points = permute(grid_points,[3 1 2]);
grid_points = grid_points(:,:);
numGridPoints = size(grid_points,2);

nProj = length(V);

x_dot_dyn_all = cell(nProj,numChoices);
x_dot_avg_all = cell(nProj,numChoices);
x_dot_combined_all = cell(nProj,numChoices);
fp_x_all = cell(nProj,numChoices);
single_trials_proj = cell(nProj,numChoices);
avg_trajectory_proj = cell(nProj,1);
unique_time_labels = unique(time_iev);
time_labels_all = [];

for iEv = 1: numel(unique_time_labels)
    
    time_labels_valid = find(time_iev == unique_time_labels(iEv));
    time_labels_all = cat(2,time_labels_all,time_labels_valid(1:end-1));
    
end


    
for iproj = 1:nProj
    
    grid_points_alignedD = V{iproj}*grid_points;
    
    for ichoice = 1: numChoices
        
        single_trials_valid =  applydynamicsSubspaceProjection(single_trials{ichoice},U');        
        avg_trajectory_valid = V{iproj}*V{iproj}'*avg_trajectory(:,:,ichoice);
        x_dot_dyn = NaN(2, numGridPoints,length(time_labels_all));
        x_dot_avg = NaN(2,numGridPoints,length(time_labels_all));
        x_dot_combined = NaN(2,numGridPoints,length(time_labels_all));
        fp_x = NaN(2,length(time_labels_all));
        assert(size(V{iproj},2) == 2, 'can only visualize a 2d plane');
        
        t_count = 1;
        
        for t_idx = 1:length(time_labels_all)
            
            resid = grid_points_alignedD - avg_trajectory_valid(:,time_labels_all(t_idx));
            prop_resid = U*A{ichoice}(:,:,time_labels_all(t_idx))*U'*resid;
            x_dot_dyn(:,:,t_count) = V{iproj}'*(prop_resid - resid);
            avg_diff = V{iproj}'*(avg_trajectory_valid(:,time_labels_all(t_idx) + 1) - avg_trajectory_valid(:,time_labels_all(t_idx)));
            
            fp_loc = inv(eye(size(U,1)) - U*A{ichoice}(:,:,time_labels_all(t_idx))*U') * ...
                (avg_trajectory_valid(:,time_labels_all(t_idx) + 1) - avg_trajectory_valid(:,time_labels_all(t_idx)));
            
            fp_x(:,t_count) = V{iproj}'*(avg_trajectory_valid(:,time_labels_all(t_idx)) + fp_loc);
            x_dot_avg(:,:,t_count) = repmat(avg_diff,[1 numGridPoints]);
            x_dot_combined(:,:,t_count) = x_dot_dyn(:,:,t_count) + x_dot_avg(:,:,t_count);
            
            t_count = t_count + 1;
            
        end
        
        x_dot_dyn_all{iproj,ichoice} = x_dot_dyn;
        x_dot_avg_all{iproj,ichoice} = x_dot_avg;
        x_dot_combined_all{iproj,ichoice} = x_dot_combined;
        fp_x_all{iproj,ichoice} = fp_x;
        avg_trajectory_proj{iproj}(:,:,ichoice) = applydynamicsSubspaceProjection(avg_trajectory_valid,V{iproj});
        single_trials_proj{iproj,ichoice} = applydynamicsSubspaceProjection(single_trials_valid,V{iproj});      
    end
    
end

analysis.V = V;
analysis.x_dot_dyn = x_dot_dyn_all;
analysis.x_dot_avg = x_dot_avg_all;
analysis.x_dot_combined = x_dot_combined_all;
analysis.single_trials = single_trials_proj;
analysis.fp_x = fp_x_all;
analysis.time_rel = time_rel(time_labels_all);
analysis.time_iev = time_iev(time_labels_all);
analysis.time = time_abs(time_labels_all);
analysis.grid_points = grid_points;
analysis.avg_trajectory_proj = avg_trajectory_proj;
end