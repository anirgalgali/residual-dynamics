function [out, out_boot] = computeDerivedDynamicsQuantities_bootstrapPairwiseAlignment(A,A_boot,sort_type,align_type,time_iev,varargin)

if(nargin == 6)
    dt = varargin{1};
end

assert(ndims(A_boot) == 4);

unique_time_labels = unique(time_iev);
V_ev = NaN(size(A));
V_ev_r = NaN(size(A));
D_ev = NaN(size(A,1),size(A,3));
U_sv = NaN(size(A));
D_sv = NaN(size(A,1),size(A,3));
V_sv = NaN(size(A));

V_ev_boot = NaN(size(A_boot));
D_ev_boot = NaN(size(A_boot,1),size(A_boot,3),size(A_boot,4));
U_sv_boot = NaN(size(A_boot));
D_sv_boot = NaN(size(A_boot,1),size(A_boot,3),size(A_boot,4));
V_sv_boot = NaN(size(A_boot));

for iEv = 1: length(unique_time_labels)

    [V_ev(:,:,time_iev == unique_time_labels(iEv)),V_ev_r(:,:,time_iev == unique_time_labels(iEv)), D_ev(:,time_iev == unique_time_labels(iEv) )]...
        = my_eigenshuffle_cmplxAdjustments(A(:,:,time_iev == unique_time_labels(iEv)), sort_type);
    
    [U_sv(:,:,time_iev == unique_time_labels(iEv)),D_sv(:,time_iev ==unique_time_labels(iEv)),...
        V_sv(:,:,time_iev ==unique_time_labels(iEv))] = my_singularshuffle(A(:,:,time_iev == unique_time_labels(iEv)),align_type); 
    
    [V_ev_boot(:,:,time_iev == unique_time_labels(iEv),:), ...
        D_ev_boot(:,time_iev == unique_time_labels(iEv),:)] = ...
        my_eigenshuffle_pairwise(D_ev(:,time_iev == unique_time_labels(iEv)),...
        V_ev(:,:,time_iev == unique_time_labels(iEv)), A_boot(:,:,time_iev == unique_time_labels(iEv),:), sort_type);
    
    [U_sv_boot(:,:,time_iev == unique_time_labels(iEv),:),...
        D_sv_boot(:,time_iev == unique_time_labels(iEv),:),...
        V_sv_boot(:,:,time_iev == unique_time_labels(iEv),:)] ...
        = my_singularshuffle_pairwise(D_sv(:,time_iev ==unique_time_labels(iEv)),...
        V_sv(:,:,time_iev ==unique_time_labels(iEv)), U_sv(:,:,time_iev == unique_time_labels(iEv)),...
        A_boot(:,:,time_iev == unique_time_labels(iEv),:), align_type);
    
end

traceCov = sum(D_sv.^2,1);
depN = sqrt(sum(D_sv.^2,1) - sum(abs(D_ev).^2,1))./sqrt(sum(abs(D_ev).^2,1));

out.A = A;
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

out_boot = repmat({out},[size(A_boot,4),1]);
for iboot = 1:size(A_boot,4)
    
    traceCov = sum(D_sv_boot(:,:,iboot).^2,1);
    depN = sqrt(sum(D_sv_boot(:,:,iboot).^2,1) - sum(abs(D_ev_boot(:,:,iboot)).^2,1))./sqrt(sum(abs(D_ev_boot(:,:,iboot)).^2,1));   
    out_boot{iboot}.A = A_boot(:,:,:,iboot);
    out_boot{iboot}.eigVal = abs(D_ev_boot(:,:,iboot));
    out_boot{iboot}.eigAngle = angle(D_ev_boot(:,:,iboot));
    [out_boot{iboot}.tau] =  compute_time_constant(out_boot{iboot}.eigVal,dt);
    [out_boot{iboot}.rotFreq] =  compute_rotation_freq(out_boot{iboot}.eigAngle,dt);
    out_boot{iboot}.eigVec = V_ev_boot(:,:,:,iboot);
    out_boot{iboot}.eigSpec = D_ev_boot(:,:,iboot);
    out_boot{iboot}.sVal = D_sv_boot(:,:,iboot);
    out_boot{iboot}.lSVec = U_sv_boot(:,:,:,iboot);
    out_boot{iboot}.rSVec = V_sv_boot(:,:,:,iboot);
    out_boot{iboot}.traceCov = traceCov;
    out_boot{iboot}.depN = depN;
    out_boot{iboot}.ev_sort_type = sort_type;
    out_boot{iboot}.sv_sort_type = align_type;
   
end



end

function [ tau] =  compute_time_constant(ev_abs,dt)
 tau = dt./log(ev_abs);
end

function [ tau] =  compute_rotation_freq(ev_ang,dt)
 tau = ev_ang./(2*pi*dt);
end