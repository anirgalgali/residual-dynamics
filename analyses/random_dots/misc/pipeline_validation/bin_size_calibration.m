%{
This script performs the bin-size calibration for the neural data to
determine the optimal bin size for residual dynamics.
Author: Aniruddh Galgali
%}
%% Set path

clearvars -except DIRS
clc

analyses_data_path = fullfile(DIRS.analysis,'/analyses/');

%% Loading Residual Data

binWidths = [15 30 45 60 90];
config_to_analyze = 3;
ref_bW = 45;

% conditionParams.targ_dir = [1;2];
opt_order = 5;
opt_lags_allBinSizes = [8;4;3;2;1];
lambda_all = 100;
animal = 'Tex';
n_dynsub_dim = 8;

%% Loading cross-validation results for the chosen monkey

analysis_name = sprintf('%s%s%d%s',animal,'_aligned_hankelagdimCV_reppasdotsTask_binsize=',ref_bW,'ms.mat');
lagdimCVresults = load(fullfile(analyses_data_path,analysis_name));


%% Fitting residual dynamics for different bin sizes

for iB = 1:length(binWidths)
       
    if(binWidths(iB) ~= ref_bW)
        analysis_name = sprintf('%s%s%d%s%d%s',animal,'_aligned_reppasdotsTask_binsize=',binWidths(iB),'ms_alignbw=',ref_bW,'ms.mat');
    else
        analysis_name = sprintf('%s%s%d%s',animal,'_aligned_reppasdotsTask_binsize=',binWidths(iB),'ms.mat');
        
    end
    load(fullfile(analyses_data_path,analysis_name));        
    nExps = length(data_aligned.residuals);

    opt_lag = opt_lags_allBinSizes(iB);
     
    for icfg = 1:nExps
     
        U_dyn = lagdimCVresults.result{icfg}.U_dyn(:,1:n_dynsub_dim);
        residual_data = data_aligned.residuals{icfg};
        trial_labels = residual_data.task_index.targ_dir;
        unique_trial_labels = unique(trial_labels);
        nconds = length(unique_trial_labels);
        
        analysis = residual_data;
        analysis = rmfield(analysis,'response');
        
        for icond = 1: nconds
            
            X = residual_data.response(:,:,trial_labels == unique_trial_labels(icond));
            X = applydynamicsSubspaceProjection(X,U_dyn);
            [A,time_idx] = estimateDynamicsTwoStageLS(X, opt_lag,lambda_all,residual_data.time_iev);
            analysis.results{icond}.A = A;
            analysis.results{icond}.time_idx = time_idx;
   
        end 
        
        analysis_all(iB,icfg) = analysis;
        
    end

end
%% Extracting rebinned eigenvalues

sort_type = 'standard';
timeWindow = {[-0.2 1.0];[-0.7 0.5]};
timeBin_refs = [15];
nConds = 2;
eigVal_ev_mean_reScaled = NaN(size(analysis.results{1}.A,1),2,length(binWidths),size(analysis_all,2),nConds);

for icfg = 1:size(analysis_all,2)
    for iB = 1:length(analysis_all)
        
        time_abs = analysis_all(iB,icfg).time;
        time_rel = analysis_all(iB,icfg).time_rel;
        time_iev = analysis_all(iB,icfg).time_iev;
        
        time_abs_dyn = time_abs(analysis_all(iB,icfg).results{1}.time_idx);
        time_rel_dyn = time_rel(analysis_all(iB,icfg).results{1}.time_idx);
        time_iev_dyn = time_iev(analysis_all(iB,icfg).results{1}.time_idx);
        
        unique_time_labels = unique(time_iev);
        
        for iEv = 1:length(unique_time_labels)
            
            idx_ev = time_iev_dyn == unique_time_labels(iEv);
            idx_valid = time_rel_dyn(idx_ev) >= timeWindow{iEv}(1) & time_rel_dyn(idx_ev) <= timeWindow{iEv}(2) ;
            
            for i_cond = 1:nConds
                
                A_ev = analysis_all(iB,icfg).results{i_cond}.A(:,:,idx_ev);
                A_ev_valid = A_ev(:,:,idx_valid);
                
                [dd_out_rescaled] =  rebin_eigenvalue_dynamics(A_ev_valid, binWidths(iB)/1000, timeBin_refs/1000, sort_type);
                eigVal_ev_mean_reScaled(:,iEv,iB,icfg,i_cond) = mean(abs(dd_out_rescaled),2);
          
                dd_out_rescaled_all{icfg,iB,iEv,i_cond} = abs(dd_out_rescaled);
            end
        end
     
    end
end


%% Setting figure properties
% 
figure;
nCols = 4;
nRows = 8;
plotWidth = 1.5;
plotHeight = 1;
figW = nCols*plotWidth;
figH = nRows*plotHeight;
fig_gain = 6;
set(gcf,'PaperUnits','centimeters')
set(gcf,'PaperSize',[figW figH])
set(gcf,'PaperPosition',[0 0 figW figH]);
set(gcf,'Units','centimeters','Position',[-44.4500   0.5556   27.1992   28.6544]);

%% Plotting rebinned eigenvalues as a function of bin size

plotpars.linewidth = 1.5;
for idim = 1:size(eigVal_ev_mean_reScaled,1)
    for icfg = 1:size(analysis_all,2)
        
        idx_plot = (idim-1)*size(analysis_all,2) + icfg;
        subplot(nRows,nCols,idx_plot);hold all;set(gca,'plotbox',[1 1 1]);
        
        for cidx = 1:2
            
            h_st_l = plot(binWidths,squeeze(eigVal_ev_mean_reScaled(idim,1,:,icfg,cidx)),'-','linewidth',plotpars.linewidth );
            h_st_m = plot(binWidths,squeeze(eigVal_ev_mean_reScaled(idim,1,:,icfg,cidx)),'o','markersize',4);
            
            h_sa_l = plot(binWidths,squeeze(eigVal_ev_mean_reScaled(idim,2,:,icfg,cidx)),'--','linewidth',plotpars.linewidth );
            h_sa_m = plot(binWidths,squeeze(eigVal_ev_mean_reScaled(idim,2,:,icfg,cidx)),'o','markersize',4);
            
            if(icfg == 1 && cidx == 1)
                
                set(h_st_l,'Color',[1 0 0]);
                set(h_st_m,'MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[1 1 1]);
                set(h_sa_l,'Color',[1 0 0]);
                set(h_sa_m,'MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[1 1 1]);
            
            elseif(icfg == 1 && cidx == 2)
                
                set(h_st_l,'Color',[0 0 1]);
                set(h_st_m,'MarkerEdgeColor',[0 0 1],'MarkerFaceColor',[1 1 1]);
                set(h_sa_l,'Color',[0 0 1]);
                set(h_sa_m,'MarkerEdgeColor',[0 0 1],'MarkerFaceColor',[1 1 1]);
                
            elseif(~(icfg == 1) && cidx == 1)
                
                set(h_st_l,'Color',[0 0 1]);
                set(h_st_m,'MarkerEdgeColor',[0 0 1],'MarkerFaceColor',[1 1 1]);
                set(h_sa_l,'Color',[0 0 1]);
                set(h_sa_m,'MarkerEdgeColor',[0 0 1],'MarkerFaceColor',[1 1 1]);
                
                
                
            elseif(~(icfg == 1) && cidx == 2)
                
                set(h_st_l,'Color',[1 0 0]);
                set(h_st_m,'MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[1 1 1]);
                set(h_sa_l,'Color',[1 0 0]);
                set(h_sa_m,'MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[1 1 1]);
                         
                
            end
                
           
        end
        
        
        if(idim == 8)
           
            set(gca,'xlim',[15 90],'xtick',[15 30 45 60 90],'xticklabel',[15 30 45 60 90]);
            xlabel(gca,'bin-size','FontName','Helvetica','FontSize',8);
            
        else
            set(gca,'xlim',[15 90],'xtick',[15 30 45 60 90],'xticklabel','');
               
        end
        
        
        if(icfg == 1)
            
            set(gca,'ylim',[0 1],'ytick',[0:0.25:1],'yticklabel',[0:0.25:1]);
            ylabel(gca,sprintf('%s%d%s','<|\lambda_{t}^{',idim,'}|>_{t}'),'FontName','Helvetica','FontSize',8);
            
        else
            set(gca,'ylim',[0 1],'ytick',[0:0.25:1],'yticklabel','');
            
        end
        
        
    end
    
end

