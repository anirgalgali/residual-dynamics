function [varargout] = plotCV_smoothness(cv_smoothness,hyper_pars,plotpars)
%{ Function that plots the cross-validated fit errors corresponding to the
% hyperparams of the second stage of the 2SLS (alpha)
% Input
% - cv_smoothness (struct)- containing the cross-validated fit error
% specific to each fold and hyperparameter setting.
% - hyper_pars (array) - specifyin all the tested hyperparam values
% - plot_pars (struct) - high-level plotting parameters
% Output
% - varargout 
%   - ah - axis handle
%   - h - figure handle
%
% Author - Aniruddh Galgali
%}

[sorted_hyper_pars,idx] = sort(hyper_pars);

switch plotpars.err_type
    
    case 'mse'
        err_train = cellfun(@(x) x.mse_train_align,cv_smoothness,'uni',false);
        err_test = cellfun(@(x) x.mse_test_align,cv_smoothness,'uni',false);

    case 'pvar'
        
        err_train = cellfun(@(x) x.pvar_train_align,cv_smoothness,'uni',false);
        err_test = cellfun(@(x) x.pvar_test_align,cv_smoothness,'uni',false);
        
end


if(ndims(err_test{1}) == 3)
   
    err_test = squeeze(cat(4,err_test{:}));
    err_test = err_test(:,:,:,idx);
    
elseif(ndims(err_test{1}) == 2)
  
    err_test = squeeze(cat(3,err_test{:}));
    err_test = err_test(:,:,idx);
    err_test = permute(err_test,[1 2 4 3]);
   
end


if(ndims(err_train{1}) == 3)
    
    err_train = squeeze(cat(4,err_train{:}));
    err_train = err_train(:,:,:,idx);
    
elseif(ndims(err_train{1}) == 2)
    
    err_train = squeeze(cat(3,err_train{:}));
    err_train = err_train(:,:,idx);
    err_train = permute(err_train,[1 2 4 3]);
    
  

end
 
mean_err_test = squeezedim(mean(mean(err_test,2),1),[1 2]);
mean_err_train = squeezedim(mean(mean(err_train,2),1),[1 2]);
std_err_train = squeezedim(std(mean(err_train,2),0,1),[1 2]);
std_err_test = squeezedim(std(mean(err_test,2),0,1),[1 2]);


switch plotpars.err_metric
    
    case 'std'
        
        err_mse_test = plotpars.err_scale * std_err_test;
        err_mse_train = plotpars.err_scale * std_err_train;
        
    case 'sem'
        
        err_mse_test = plotpars.err_scale * (1/sqrt(size(std_err_test,2))).*std_err_test;
        err_mse_train = plotpars.err_scale * (1/sqrt(size(std_err_test,2))).*std_err_train;
        
end

figure;
nalign = size(err_test,3);
for  ialign = 1:nalign
   
    ah = subplot(1,nalign,ialign);hold(ah);set(ah,'plotbox',[1 1 1]);
    hh_test = ploterr(log10(sorted_hyper_pars), mean_err_test(ialign,:),[],err_mse_test(ialign,:),'-','abshhy', plotpars.bar_length);
    set(hh_test(1),'Color',plotpars.col_test,'linewidth',plotpars.linewidth);
    set(hh_test(2),'Color',plotpars.col_test,'linewidth',plotpars.bar_width);
    set(hh_test,'Parent',ah);
    hm_test = plot(ah,log10(sorted_hyper_pars), mean_err_test(ialign,:),plotpars.marker_type);
    set(hm_test,'markeredgecolor',plotpars.col_test ,'markersize',plotpars.markersize,'markerfacecolor',[1 1 1]);
    
    if(plotpars.do_plot_train)
        
        hh_train = ploterr(log10(sorted_hyper_pars), mean_err_train(ialign,:),[],err_mse_train(ialign,:),'-','abshhy', plotpars.bar_length);
        set(hh_train(1),'Color',plotpars.col_train,'linewidth',plotpars.linewidth);
        set(hh_train(2),'Color',plotpars.col_train,'linewidth',plotpars.bar_width);
        
        set(hh_train,'Parent',ah);
        
        hm_train = plot(ah,log10(sorted_hyper_pars), mean_err_train(ialign,:),plotpars.marker_type);
        set(hm_train,'markeredgecolor',plotpars.col_train ,'markersize',plotpars.markersize,'markerfacecolor',[1 1 1]);
        
    end
    
    xlabel('log_{10}(\alpha)')
    if(ialign == 1)
        ylabel(sprintf('%s%s','cv error ',plotpars.err_type))
    end
end

h = gcf;
varargout{1} = h;

end