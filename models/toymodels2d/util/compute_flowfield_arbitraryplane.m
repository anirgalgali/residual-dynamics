function [locs, flow] = compute_flowfield_arbitraryplane(A, U, V, gridX, gridY, dynamics_type, varargin)
%{ This function computes a flow-field in ar arbitrary state-space plane.
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
% - dynamics_type (string) - specifies whether the dynamics specified in A
%   is discrete-time or continuous time. 
% - varargin (float) - indicate the step-size (dt) between consecutive time
%   points in A
%
% Output
% - analysis (struct) - containing the estimates of the flow-field.
%} 

if(strcmp(dynamics_type,'discrete_time'))
    assert(nargin == 7,'need step size argument')
    step_size = varargin{1};
end
[grid_locs_X,grid_locs_Y] = meshgrid(gridX,gridY);

locs_ = cat(3,grid_locs_X,grid_locs_Y);
locs_ = permute(locs_,[3 1 2]);

grid_points_all = locs_(:,:);
locs = V * grid_points_all;
locs = reshape(V'*locs,[size(V,2), size(locs_,2),size(locs_,3)]);


switch dynamics_type
    case 'continuous_time'
        flow = V' * U * A * U'* V * grid_points_all;
    case 'discrete_time'
        flow = V' *(U * A * U' - eye(size(U,1)))* V * grid_points_all;
        flow = (1./step_size)*flow;
        
end
flow = reshape(flow,[size(V,2), size(locs,2), size(locs,3)]);


end
