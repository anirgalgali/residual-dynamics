function [out] = computeDerivedDynamicsQuantities(A,sort_type,align_type,time_iev,dt,varargin)
%{ This function computes the properties of time-varying dynamics matrices 
%  that are estiamted using the residual dynamics pipeline.
% Inputs
% - A (3d tensor) - of size n_dyn_dim x n_dyn_dim x n_times, where n_dyn_dim
%  = dimensionality of a sequence (n_times) of dynamics matrices
% - sort_type (string) - that determines how to sort eigenvalues at the first
%   time step (typically 'abs' - meaning that at the first time eigenvalues
%   are sorted according to their absolute values)
% - align_type (string) - determines how to sort the singular values -
%   whether based on the right or left singular vectors.
% - time_iev - time-labels associated with each time-bin of A, which indicate
%   which temporal alignment a given A matrix belongs to
% - dt - time-step (i.e difference in times between two adjacent A matrices)
% -varargin
%  - order_across_alignment (boolean) - true/false indicates whether to
%    retain ordering of time-varying eigen/singular values across temporal
%    alignments.
%
% Outputs
% -out (struct) - containing the various properties of the residual
%                dynamics defined by A
%
%
%
% Author : Aniruddh Galgali (Written : Apr 2019)
%}
if(nargin == 5)
    order_across_alignment = false;
elseif(nargin == 6)  
    order_across_alignment = varargin{1};
end

unique_time_labels = unique(time_iev);
V_ev = NaN(size(A));
V_ev_r = NaN(size(A));
D_ev = NaN(size(A,1),size(A,3));
U_sv = NaN(size(A));
D_sv = NaN(size(A,1),size(A,3));
V_sv = NaN(size(A));

if(~order_across_alignment)
    for iEv = 1: length(unique_time_labels)

        [V_ev(:,:,time_iev == unique_time_labels(iEv)),V_ev_r(:,:,time_iev == unique_time_labels(iEv)), D_ev(:,time_iev == unique_time_labels(iEv) )]...
            = my_eigenshuffle(A(:,:,time_iev == unique_time_labels(iEv)), sort_type);

        [U_sv(:,:,time_iev == unique_time_labels(iEv)),D_sv(:,time_iev ==unique_time_labels(iEv)),...
            V_sv(:,:,time_iev ==unique_time_labels(iEv))] = my_singularshuffle(A(:,:,time_iev == unique_time_labels(iEv)),align_type); 

    end
    
else
    
    [V_ev,V_ev_r,D_ev] = my_eigenshuffle(A, sort_type);
    [U_sv,D_sv,V_sv] = my_singularshuffle(A,align_type);

    
end

traceCov = sum(D_sv.^2,1);
depN = sqrt(sum(D_sv.^2,1) - sum(abs(D_ev).^2,1))./sqrt(sum(abs(D_ev).^2,1));

out.eigVal = abs(D_ev);
out.eigAngle = angle(D_ev);
[out.tau] =  compute_time_constant(out.eigVal,dt);
[out.rotFreq] =  compute_rotation_freq(out.eigAngle,dt);
out.eigVec = V_ev;
out.eigVec_r = V_ev_r;
out.eigSpec = D_ev;
out.sVal = D_sv;
out.lSVec = U_sv;
out.rSVec = V_sv;
out.traceCov = traceCov;
out.depN = depN;
out.ev_sort_type = sort_type;
out.sv_sort_type = align_type;

end

function [ tau] =  compute_time_constant(ev_abs,dt)
 tau = dt./log(ev_abs);
end

function [ tau] =  compute_rotation_freq(ev_ang,dt)
 tau = ev_ang./(2*pi*dt);
end

