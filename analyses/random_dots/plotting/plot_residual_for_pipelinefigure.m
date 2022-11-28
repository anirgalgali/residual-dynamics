clearvars -except DIRS
clc
rseed = rng('default');
close all
%%

pexp_residuals = ProjectLoad('array_dotsTask_alignedResiduals_dim=20_bW=45ms_doFR=0_doSqrt=1_nonSmooth');
pexp_lagdimCV = ProjectLoad_v2('array_dotsTask_LagDimHankelCV_dim=20_bW=45ms_doFR=0_doSqrt=1_nonSmooth');
pexp_singleTrials = ProjectLoad_v2('array_dotsTask_alignedTrialsPooled_dim=20_bW=45ms_doFR=0_doSqrt=1_nonSmooth');
animal = 'Tex';
%%
opts_final.do_cross_validate_noise_model = false;
n_dim = 8;
lag = 3;
alpha_tex = [200;50];
alpha_vito = alpha_tex;
% alpha_vito = [5000;500];
opts_final.dim = {n_dim;n_dim};
opts_final.lag = {lag;lag};
opts_final.order_across_alignment = true;
% result_final = cell(1,n_configs);
cidx_to_plot = 1;
switch animal 
    case 'Tex'
        opts_final.alpha = {[alpha_tex];[alpha_tex]};
        
    case 'Vito'
        opts_final.alpha = {[alpha_vito];[alpha_vito]};
end
id_list = pexp_residuals.idlist;
idx_animal_valid = cell2mat(cellfun(@(x) ~isempty(x),cellfun(@(x) strfind(x,animal),id_list,'uni',false),'uni',false));
id_list_selected = id_list(idx_animal_valid);
config_idx = 3;
trial_idx = 2000;
ev_idx = 1;
for iexp = config_idx
    
    result_lagdim_cv = AnalysisLoad_v2(pexp_lagdimCV.name, id_list_selected{iexp});
    analysis = AnalysisLoad(pexp_residuals.name,id_list_selected{iexp});
    time_labels = analysis.residuals.time_iev;
    trial_labels = analysis.residuals.task_index.targ_dir;
    data_residuals = analysis.residuals.response;
    opts_final.hankel_order = result_lagdim_cv.fit_pars.hankel_order;
    opts_final.hankel_rank = result_lagdim_cv.opt_hankel_rank;
    
    [result_final] = fit_final_resdyn_model(data_residuals, ...
        analysis.residuals.time, analysis.residuals.time_rel,...
        time_labels, trial_labels, opts_final);
    idx_ev  = time_labels == ev_idx;       
    [data_fs_in] = applydynamicsSubspaceProjection(data_residuals(:,idx_ev,trial_labels == cidx_to_plot),...
        result_final.U_dyn(:,1:n_dim));
    
    figure;set(gcf,'position', [90  26  498 1068])

    for idim = 1:n_dim
        ah = subplot(n_dim,1,idim);hold(ah)
        plot(squeeze(data_fs_in(idim,lag+1:end,trial_idx)));
        set(ah,'ylim',[-2 2],'ytick',[-2 0 2]);

 
    end
    
    figure;set(gcf,'position', [90  26  498  1068])

    for idim = 1:n_dim
        ah = subplot(n_dim,1,idim);hold(ah)
        plot(squeeze(result_final.final_model{cidx_to_plot}.X_fs{ev_idx}(idim,:,trial_idx)));
        set(ah,'ylim',[-1 1],'ytick',[-1 0 1]);

    end
    
end
