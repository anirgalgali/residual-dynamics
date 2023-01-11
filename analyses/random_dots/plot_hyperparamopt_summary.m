%{This script plots the summary of the optimal hyperparameters across both 
% monkeys for each task configuration
% Author: Aniruddh Galgali (Apr 2022)
%}

clearvars -except DIRS
clc
rseed = rng('default');

%% Choose animal and bin-size
animal = {'Tex','Vito'};
bin_size = 45;
data_path = './data/analyses/';
result_cvlagdim = cell(length(animal),1);
result_cvsmooth = cell(length(animal),1);

for iani = 1:length(animal)
   file_name = sprintf('%s%s%d%s',animal{iani},'_aligned_hankelagdimCV_reppasdotsTask_binsize=',bin_size,'ms.mat');
   R = load(fullfile(data_path,file_name));
   result_cvlagdim{iani} = R.result;
   
   file_name = sprintf('%s%s%d%s',animal{iani},'_aligned_smoothnessCV_reppasdotsTask_binsize=',bin_size,'ms.mat');
   R = load(fullfile(data_path,file_name));
   result_cvsmooth{iani} = R.result_smoothness;
   
end
%%
col_mat = [0 0 1; 1 0 0];
marker_style = {'^';'d';'o';'s'};
marker_size = 10;
%%
figure;
for iani = 1: length(animal)
    cv_lagdim_opt = cell(length(result_cvlagdim{iani}),1);
    dims_all = cell(1, length(result_cvlagdim{iani}));
    lags_all = cell(1, length(result_cvlagdim{iani}));
    for iexp = 1:length(result_cvlagdim{iani})
        
        cv_lagdim_opt{iexp} = result_cvlagdim{iani}{iexp}.cv_lagdim_opt;
        dims_all{iexp} = cell2mat(cellfun(@(y) y.mse_opt.dims(1), cv_lagdim_opt{iexp},'uni',false));
        lags_all{iexp} = cell2mat(cellfun(@(y) y.mse_opt.lags(1), cv_lagdim_opt{iexp},'uni',false));
    end
    dims_all = cell2mat(dims_all);
    lags_all = cell2mat(lags_all);

    ah = subplot(1,length(animal),iani);set(ah,'plotbox',[1 1 1]);hold all;
    set(gca,'ylim',[1.5 12.5],'ytick',[2:2:12],'yticklabel',[2:2:12],'xlim',[0.5 2.5],'xtick',[1 2],'xticklabel',{'dim';'lag'});
    
    for iexp = 1:size(dims_all,2)
        
        if(iexp == 1)
            dims_all(:,iexp) = flipud(dims_all(:,iexp));
            lags_all(:,iexp) = flipud(lags_all(:,iexp));
        end
        
        for icond = 1:size(dims_all,1)
            plot(1+ 0.04*randn(1),dims_all(icond,iexp) + 0.04*randn(1),marker_style{iexp},'markersize',marker_size,'MarkerFaceColor',col_mat(icond,:),'MarkerEdgeColor',[ 0 0 0]);
            plot(2+ 0.04*randn(1),lags_all(icond,iexp) + 0.04*randn(1),marker_style{iexp},'markersize',marker_size,'MarkerFaceColor',col_mat(icond,:),'MarkerEdgeColor',[ 0 0 0]);
        end
    end
    title(sprintf('%s%s',animal{iani},': optimal lag and dim (1 s.e.m rule)'))
end

%%
figure;
for iani = 1: length(animal)
    cv_smooth_opt = cell(length(result_cvsmooth{iani}),1);
    alpha_all_max = cell(1, length(result_cvsmooth{iani}));
    alpha_all_min = cell(1, length(result_cvsmooth{iani}));
    for iexp = 1:length(result_cvlagdim{iani})
        cv_smooth_opt{iexp} = result_cvsmooth{iani}{iexp}.cv_smooth_opt;
        alpha_all_max{iexp} = cell2mat(cellfun(@(y) y.T_mse.alpha_mse(1),cv_smooth_opt{iexp},'uni',false));
        alpha_all_min{iexp} = cell2mat(cellfun(@(y) y.T_mse.alpha_mse(end),cv_smooth_opt{iexp},'uni',false));
    end

    alpha_all_max = cat(3,alpha_all_max{:});
    alpha_all_min = cat(3,alpha_all_min{:});
    
    ah = subplot(1,length(animal),iani);set(ah,'plotbox',[1 1 1]);hold all;
    set(gca,'ylim',[-1.5 5.5],'ytick',[-1:1:5],'yticklabel',[-1:1:5],'xlim',[0.5 2.5],'xtick',[1 2],'xticklabel',{'stim';'sacc'});
    ylabel(gca,'log_{10}(\alpha)')
    for iexp = 1:size(alpha_all_max,3)
        
        if(iexp == 1)
            alpha_all_max(:,:,iexp) = fliplr(alpha_all_max(:,:,iexp));
            alpha_all_min(:,:,iexp) = fliplr(alpha_all_min(:,:,iexp));
        end
        
        for icond = 1:size(dims_all,1)
            plot(1+ 0.04*randn(1),log10(squeeze(alpha_all_max(1,icond,iexp))) + 0.02*randn(1),marker_style{iexp},'markersize',marker_size,'MarkerFaceColor',col_mat(icond,:),'MarkerEdgeColor',col_mat(icond,:));
            plot(2+ 0.04*randn(1),log10(squeeze(alpha_all_max(2,icond,iexp))) + 0.02*randn(1),marker_style{iexp},'markersize',marker_size,'MarkerFaceColor',col_mat(icond,:),'MarkerEdgeColor',col_mat(icond,:));
            
            plot(1+ 0.04*randn(1),log10(squeeze(alpha_all_min(1,icond,iexp))) + 0.02*randn(1),marker_style{iexp},'markersize',marker_size,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',col_mat(icond,:));
            plot(2+ 0.04*randn(1),log10(squeeze(alpha_all_min(2,icond,iexp))) + 0.02*randn(1),marker_style{iexp},'markersize',marker_size,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',col_mat(icond,:));
        end
    end
    title(sprintf('%s%s',animal{iani},': optimal smoothness reg.const'))

end
