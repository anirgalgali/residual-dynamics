function [varargout] = plotBehaviour(animal,data_path,plotpars)
% Loading data

pexp_beh = ProjectLoad_v2(sprintf('%s%s','array_dotsTask_beh_',animal));
load(fullfile(data_path, sprintf('%s%s',animal,'_aligned_reppasdotsTask_binsize=45ms.mat')));

idlist_align_animal = [data_aligned.expt_names{:}];
idlist_align_animal = cellfun(@(x) strrep(x,'-BhvDspEFit.mat',''),idlist_align_animal,'uni',false);
valid_sessions = false(length(pexp_beh.idlist),1);
for ii = 1 : length(pexp_beh.idlist)
    
    analysis{ii} = AnalysisLoad_v2(pexp_beh.name,pexp_beh.idlist{ii},'fast');
    
    if(~isfield(analysis{ii}.targInfo,'config'))
        
        allTargetConfigs(ii) = NaN;
    
    else
        allTargetConfigs(ii) = analysis{ii}.targInfo.config;
        valid_sessions(ii) = contains(pexp_beh.idlist{ii},idlist_align_animal);
        
    end
     
    analysis{ii}.behaviour.percChoiceToPref;
    
end

% valid_sessions = ~isnan(allTargetConfigs) & ismember(pexp_beh.idlist,idlist_align_animal);
invalid_file_names = cellfun(@(x) x.thisFile,analysis(~valid_sessions),'uni',false);
analysis = analysis(valid_sessions);

plotPsychometricData(analysis,plotpars)
set(gcf,'Position',[20    63   690   576]);

if(plotpars.doSave)
    if(plotpars.plotAllData)
        figName = sprintf('%s%s%s%s','psychometricData_',animal,'_allSegments_','separateConfigSigmoidFits');
    else
        figName = sprintf('%s%s%s','psychometricData_',animal,'separateConfigSigmoidFits');
    end
    export_fig(fullfile(plotpars.save_fig_path,figName),'-painters','-transparent','-pdf')
end

plotTargetConfigurations(analysis,plotpars)
set(gcf,'Position',[-5  408 1265  259]);
    
if(plotpars.doSave)

    figName = sprintf('%s%s%s','configOnScreen_',animal,'_allSessions');
    export_fig(fullfile(plotpars.save_fig_path,figName),'-painters','-transparent','-pdf')
end

varargout{1} = invalid_file_names;

end


function plotPsychometricData(analysis,plotpars)

numTrials_allSegments = [];
targetConfigs_allSegments = [];
percChoiceToPref_allSegments = {};
dotsCoh_allSegments = {};
n_count = 1;

for ii = 1: length(analysis)
    numSegments = length(analysis{ii}.segments.valid);
    for kk = 1: numSegments
        if(analysis{ii}.segments.valid(kk))
            
            numTrials_allSegments(n_count) = sum(analysis{ii}.segments.idxs{kk});
            targetConfigs_allSegments(n_count) = analysis{ii}.targInfo.config;
            percChoiceToPref_allSegments(n_count) = analysis{ii}.behaviour.percChoiceToPref(kk);
            dotsCoh_allSegments(n_count) = analysis{ii}.behaviour.signed_coh(kk);
            n_count = n_count + 1;
            
        end
    end
end

valid_segments = ~isnan(targetConfigs_allSegments);
targetConfigs_allSegments = targetConfigs_allSegments(valid_segments);
numTrials_allSegments = numTrials_allSegments(valid_segments);
[~,ind] = sort(numTrials_allSegments,'descend');
targetConfigs_allSegments_sorted  = targetConfigs_allSegments(ind);

percChoiceToPref_allSegments = percChoiceToPref_allSegments(valid_segments);
dotsCoh_allSegments = dotsCoh_allSegments(valid_segments);

[uniqueTargetConfigs,~,idx] = unique(targetConfigs_allSegments );
marker_cols_all = NaN(length(percChoiceToPref_allSegments),3);
marker_type_all = cell(length(percChoiceToPref_allSegments),1);
for jj = 1: length(uniqueTargetConfigs)
    marker_cols_all(idx == jj,:) = repmat(plotpars.markercol(jj,:),[sum(idx == jj) 1]);
    largest_segment_idx(jj) = ind(find(targetConfigs_allSegments_sorted == uniqueTargetConfigs(jj),1));
    marker_type_all(idx == jj) = repmat(plotpars.marker_type(jj),[sum(idx == jj) 1]) ;
end

%% Plotting all the data across all sessions
figure; ah = gca; hold(ah);set(ah,'plotbox',[1 1 1]);

if(plotpars.plotAllData)
    for ii = 1: length(percChoiceToPref_allSegments)
        
        
        
        xx = dotsCoh_allSegments{ii};
        yy = percChoiceToPref_allSegments{ii};
        
        hh = scatter(xx,yy,plotpars.markersize,marker_type_all{ii});
        
        set(hh,'markeredgecolor',marker_cols_all(ii,:),'markerfacecolor'...
            ,marker_cols_all(ii,:),'markerFaceAlpha',plotpars.marker_face_alpha,'markerEdgeAlpha',plotpars.marker_edge_alpha);
        
        set(gca,'plotbox',[1 1 1]);
        
        
    end
    
end

%% Plotting line
switch plotpars.line_type
    
    case 'segment_average'
        
        mean_percPref = cell(length(uniqueTargetConfigs),1) ;
        
        for jj = 1: length(uniqueTargetConfigs)
            
            idx_config = targetConfigs_allSegments == uniqueTargetConfigs(jj);
            all_dots_coh = cell2mat(dotsCoh_allSegments(idx_config)');
            all_percPref = cell2mat(percChoiceToPref_allSegments(idx_config)');
            [all_unique_dots_coh] = unique(all_dots_coh);
        
            for ii = 1: length(all_unique_dots_coh)
                
                idx_avg = all_dots_coh == all_unique_dots_coh(ii) ;
                mean_percPref{jj}(ii) = mean(all_percPref(idx_avg));
                
            end
            
            mean_percPref{jj} = smooth(mean_percPref{jj},3);
            hh_l = plot(all_unique_dots_coh,mean_percPref{jj},'-');
            hh_m = plot(all_unique_dots_coh,mean_percPref{jj},plotpars.marker_type{jj},'markersize',8);
            set(hh_l,'Color',plotpars.markercol(jj,:),'linewidth',4)
            set(hh_m,'markerfacecolor',[1 1 1],'markeredgecolor',plotpars.markercol(jj,:));
        
        end
       
        
        
    case 'fit_logistic'
        
        for jj = 1: length(uniqueTargetConfigs)
            
            idx_config = targetConfigs_allSegments == uniqueTargetConfigs(jj);
            all_dots_coh = cell2mat(dotsCoh_allSegments(idx_config)');
            all_percPref = cell2mat(percChoiceToPref_allSegments(idx_config)');
            
            
            ft = fittype( 'logisticFunc( x, k, x0 )' );
            ff = fit( all_dots_coh, all_percPref, ft, 'StartPoint', [2,0] );
            hh_l = plot(ff,'-');
            set(hh_l,'Color',plotpars.markercol(jj,:),'linewidth',2)
            hh_m = plot(hh_l.XData(1:50:end),hh_l.YData(1:50:end),plotpars.marker_type{jj});
            set(hh_m,'markerfacecolor',[1 1 1],'markeredgecolor',plotpars.markercol(jj,:),'markersize',8);
            
        end
        
        if(plotpars.plotBinnedAverage)
            coh_bins = plotpars.coh_bins ;
            mean_percPref = cell(length(uniqueTargetConfigs),1) ;
            
            for jj = 1: length(uniqueTargetConfigs)
                
                
                idx_config = targetConfigs_allSegments == uniqueTargetConfigs(jj);
                all_dots_coh = cell2mat(dotsCoh_allSegments(idx_config)');
                all_percPref = cell2mat(percChoiceToPref_allSegments(idx_config)');
                [N,edges,bin_idx] = histcounts(all_dots_coh,coh_bins);
                [all_unique_bin_idx] = unique(bin_idx);
                
                for ii = 1: length(all_unique_bin_idx)
                    
                    idx_avg = bin_idx == all_unique_bin_idx(ii) ;
                    mean_percPref{jj}(ii) = mean(all_percPref(idx_avg));
                    std_percPref{jj}(ii) = std(all_percPref(idx_avg));
                    sem_percPref{jj}(ii) = (1./sqrt(sum(idx_avg))).*std(all_percPref(idx_avg));
                end
                
                bin_centers = 0.5*(coh_bins(1:end-1) + coh_bins(2:end));
                hh_m = errorbar(bin_centers(all_unique_bin_idx)',mean_percPref{jj}', 2.*sem_percPref{jj}',plotpars.marker_type{jj});
                set(hh_m,'Color',plotpars.markercol(jj,:),'markeredgecolor',plotpars.markercol(jj,:),...
                    'markerfacecolor',[1 1 1],'markersize',10,'linewidth',2);
            end
            
        end
        
end

set(gca,'plotbox',[1 1 1]);
axis tight
xlabel('signed coherency (%)')
ylabel('fraction of choices to pref (a.u)')
set(gca,'FontName',plotpars.FontName,'FontSize',plotpars.FontSize)
legend off

end


function plotTargetConfigurations(analysis,plotpars)

[allTargConfigs] = cellfun(@(x) x.targInfo.config,analysis,'uni',false);

allTargConfigs = cell2mat(allTargConfigs);
valid_sessions =~isnan(allTargConfigs);
allTargConfigs = allTargConfigs(valid_sessions);
unique_targ_configs = unique(allTargConfigs);
figure;
% pf = bda_figure('',[length(unique_targ_configs) 1],1);
analysis = analysis(valid_sessions);
for jj = 1:length(unique_targ_configs)
    
    ah = subplot(1,length(unique_targ_configs),jj);
    idx = allTargConfigs == unique_targ_configs(jj);
    analysis_valid = analysis(idx);
    hold(ah);set(ah,'plotbox',[1 1 1]);
    plot(0,0,plotpars.marker_type{jj},'markersize',plotpars.fp_markersize,'Color',plotpars.fp_col(jj,:));
    
    for ii = 1 : length(analysis_valid)
        
        h_T1 = plot(analysis_valid{ii}.targInfo.relPos_T1(1) + 0.5*randn(1), analysis_valid{ii}.targInfo.relPos_T1(2) + 0.5*randn(1),'o');
        h_T2 = plot(analysis_valid{ii}.targInfo.relPos_T2(1) + 0.5*randn(1), analysis_valid{ii}.targInfo.relPos_T2(2) + 0.5*randn(1),'o');
        
        
        if(strcmp(analysis_valid{ii}.targInfo.labels(analysis_valid{ii}.targInfo.pref_idx),'T1'))
           
            set(h_T1,'markersize',plotpars.targ_markersize,'MarkerFaceColor',plotpars.pref_col,'MarkerEdgeColor',plotpars.pref_col)
            set(h_T2,'markersize',plotpars.targ_markersize,'MarkerFaceColor',plotpars.anti_col,'MarkerEdgeColor',plotpars.anti_col)
            
            
        elseif(strcmp(analysis_valid{ii}.targInfo.labels(analysis_valid{ii}.targInfo.pref_idx),'T2'))
            
            set(h_T2,'markersize',plotpars.targ_markersize,'MarkerFaceColor',plotpars.pref_col,'MarkerEdgeColor',plotpars.pref_col)
            set(h_T1,'markersize',plotpars.targ_markersize,'MarkerFaceColor',plotpars.anti_col,'MarkerEdgeColor',plotpars.anti_col)
            
            
        end
        
        
        
        
    end
    
    
    set(gca,'plotbox',[1 1 1])
    set(gca,'xlim',[-15 15],'ylim',[-15 15],'xtick',[0],'ytick',0,'xticklabel',[0],'yticklabel',0)
    xlabel('screen X')
    ylabel('screen Y')
    set(gca,'FontName',plotpars.FontName,'FontSize',plotpars.FontSize)
    
end


end

