function [out] = compute_eigenvector_taskrelaxes_overlap(evec, eval, U_dyn, task_rel_axes, dim_masks, vars_to_store)
   
% Computes overlap between task relevant vectors of the two-area model
% and the estimated residual eigenvectors at a given time 't'.

%INPUTS
% eval - n_dim x 1 vector of eigenvalues (n_dim = dim of dynamics)
%
% evec - n_dim x n_dim matrix, with columns being the eigenvectors
%
% U_dyn - obs_dim x obs_dim matrix, with columns being the dimensions
%         that define the dynamics subspace
%
% task_rel_axes - nested struct containing task relevatn model directions 
%                  organized by area and direction type. First field
%                  level specifies area - {'ppc' or 'pfc'} and second
%                  field level specifies direction type - {'choice' or 'time'}
%
% dim_masks - obs_dim x 1 logical vector that specifies how dimensions of 
%             the task relevant axes get partioned into areas 
%
% vars_to_store - cell array specifying output variables to store

%OUTPUTS
% out - nested struct with same organization as variable task_rel_axes,
%       with the output variables specified in variable vars_to_store 

% Author: Aniruddh Galgali, Apr 2021


%Extract names of individual modules that make up the network
area_names = fieldnames(task_rel_axes);


unique_eval = unique(eval);
n_dim = length(eval);
evec = U_dyn(:,1:n_dim)*evec;

for iev = 1:length(unique_eval)

    idx_valid = eval == unique_eval(iev);

    for iarea = 1:length(area_names)

        direction_names = fieldnames(task_rel_axes.(area_names{iarea}));

        for idir = 1: length(direction_names)
        
            % Extracting the task-relevatn vector 
               
            task_rel_axis = task_rel_axes.(area_names{iarea}).(direction_names{idir})(dim_masks);
            
            % Computing projection of task-relevant vector with each unique
            % eigenvector
            if(norm(task_rel_axis) ~= 0)
                proj_val =  normc(task_rel_axis)'*  normc(real(evec(:,idx_valid)));
            else
                proj_val =  task_rel_axis'*  normc(real(evec(:,idx_valid))); % if task_rel_axis is a zero vector, do not normalize before computing projection
            end   

            for ivar = 1:length(vars_to_store)

                switch vars_to_store{ivar}
                    
                    % Storing raw projection
                    case 'proj'

                        out.(area_names{iarea}).(direction_names{idir}).(vars_to_store{ivar})(idx_valid) = proj_val;
                    
                    % Storing angle (degrees) corresponding to projection
                    case 'angle'

                        angle_val = rad2deg(acos(proj_val));
                        if(any(angle_val > 90))
                            idx_inv = angle_val > 90;
                            angle_val(idx_inv) = 180 - angle_val(idx_inv);
                        end
                        out.(area_names{iarea}).(direction_names{idir}).(vars_to_store{ivar})(idx_valid) = angle_val;

                end

            end

        end
    end
end
end