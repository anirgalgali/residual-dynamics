function [data_agg] =  compute_ovelap_resdyn_avgtaskplanes(analysis_condavg, result_final)
%{ This function computes the overlap between the eigenvectors of the
%  residual dynamics and the task planes defined using aligned,
%  trial-averaged responses. The overlap is computed separately for each
%  task configuration and then aggregated across configurations.
%  Input
%  -analysis_condavg (cell array) - of size n_configs x 1, each cell 
%   contains the aligned task planes and the associated information
%  -result_final - (cell array) - of size n_configs x 1, each cell 
%   contains the residual dynamics fits (of all choice conditions/task epochs
%   of the corresponding configuration
%  Outputs
%  - data_agg (table) - indicaitng the different ev/sv values, the overlap
%  of the corresponding eigenvectors and task planes,...
%
%  Author - Aniruddh Galgali (May 2021)
%}
data_agg = struct();
n_configs = length(analysis_condavg);
assert(length(result_final) == n_configs,'mismathc in number of data points')

c_count = 1;

for iexp = 1:n_configs
    
    nconds = length(result_final{iexp}.final_model);
    alignments = fieldnames(analysis_condavg{iexp});
    direction_labels = analysis_condavg{iexp}.(alignments{1}).direction_labels;
    
    for icond = 1:nconds
        
        for tt = 1:size(result_final{iexp}.final_model{icond}.eigVal,2)
            
            eig_val = result_final{iexp}.final_model{icond}.eigVal(:,tt);
            [max_eig_val, max_eig_val_idx] = max(eig_val);
            
            eig_vec_t_real = real(result_final{iexp}.final_model{icond}.eigVec_r(:,:,tt));            
            eig_angle = result_final{iexp}.final_model{icond}.eigAngle(:,tt);
            rot_freq = result_final{iexp}.final_model{icond}.rotFreq(:,tt);     
            idx_imag = eig_angle ~= 0;
            
            eig_vec_t_full_real = result_final{iexp}.U_dyn(:,1:size(eig_vec_t_real,1)) * eig_vec_t_real;
            
            time_label_t = result_final{iexp}.final_model{icond}.time_iev(tt);
            
            unique_abs_eig_val = unique(eig_val);
            
            for id = 1: length(unique_abs_eig_val) 
                
                idx_ev = eig_val == unique_abs_eig_val(id); 
                
                for idir = 1: length(direction_labels)
                    
                    idx_dir = ismember(analysis_condavg{iexp}.(alignments{time_label_t}).direction_labels,...
                        direction_labels{idir});
                      
                    all_angles.(direction_labels{idir})(idx_ev) = ...
                                rad2deg(subspacea(analysis_condavg{iexp}.(alignments{time_label_t}).directions{idx_dir},...
                                eig_vec_t_full_real(:,idx_ev)));
                            
                  
                  
                end
                
            end
            
   
            for id = 1:size(eig_val,1)
                
                data_agg(c_count).ev = eig_val(id);
                data_agg(c_count).ev_idx = id;
                data_agg(c_count).max_ev = max_eig_val;
                data_agg(c_count).rot_freq = abs(rot_freq(id));
                
                if(id == max_eig_val_idx)
                    data_agg(c_count).is_max = 1;
                else
                    data_agg(c_count).is_max = 0;
                end
                
                if(idx_imag(id))
                    data_agg(c_count).is_imag = 1;
                else
                    data_agg(c_count).is_imag = 0;
                end
                
                plane_names = fieldnames(all_angles);
                for iplane = 1: length(plane_names)
                    data_agg(c_count).(['angle_with_' plane_names{iplane}]) = ...
                        all_angles.(plane_names{iplane})(id);
                    
                end
 
                data_agg(c_count).config_idx = iexp;
                data_agg(c_count).time = result_final{iexp}.final_model{icond}.time(tt);
                data_agg(c_count).time_rel = result_final{iexp}.final_model{icond}.time_rel(tt);
                data_agg(c_count).time_label = result_final{iexp}.final_model{icond}.time_iev(tt);
                data_agg(c_count).time_idx = tt;
                data_agg(c_count).choice_idx = icond;
                
                c_count = c_count + 1;
                
                
            end
            
            
        end
        
    end
    
end
data_agg = struct2table(data_agg);

end