function [varargout] = plotCV_lagDim(cv_joint,hyper_pars,plotpars)
%{ Function that plots the cross-validated fit errors corresponding to the
% hyperparams of the first stage of the 2SLS (lag and dime)
% Input
% - cv_joint(struct)- containing the cross-validated fit error
% specific to each fold and hyperparameter setting.
% - hyper_pars (array) - specifying all the tested hyperparam values
% - plot_pars (struct) - high-level plotting parameters
% Output
% - varargout 
%   - ah - axis handle
%   - h - figure handle
%
% Author - Aniruddh Galgali
%}

assert(istable(hyper_pars),'hyper_pars should be of type table');

switch plotpars.err_type
    
    case 'mse'
        
        err_test = cellfun(@(x) x.mse_test,cv_joint,'uni',false);
        err_train = cellfun(@(x) x.mse_train,cv_joint,'uni',false);
        
    case 'pvar'
        err_test = cellfun(@(x) x.pvar_test,cv_joint,'uni',false);
        err_train = cellfun(@(x) x.pvar_train,cv_joint,'uni',false);
end

[n_cv,n_folds,n_conds] = size(err_test{1});

unique_dims = unique(hyper_pars.dim);
unique_lags = unique(hyper_pars.lag);

if(n_cv > 1)
    idx_cv_plot = randperm(n_cv);
    fprintf('multipe cv runs present. Choosing run %d to plot\n',idx_plot)
else
    idx_cv_plot = 1;
end

figure;

for i_cond = 1:n_conds
    
    ah = subplot(n_conds,1,i_cond);hold(ah);set(ah,'plotbox',[1 1 1]);
    
    err_test_cond_mean = cell2mat(cellfun(@(x) mean(x(:,:,i_cond),2)', err_test,'uni',false));
    err_test_cond_std = cell2mat(cellfun(@(x) std(x(:,:,i_cond),0,2)', err_test,'uni',false));
    err_train_cond_mean = cell2mat(cellfun(@(x) mean(x(:,:,i_cond),2)', err_train,'uni',false));
    err_train_cond_std = cell2mat(cellfun(@(x) std(x(:,:,i_cond),0,2)', err_train,'uni',false));
    
    switch plotpars.x_type
        
       
        case 'lag'
            
            assert(size(plotpars.col_test,1) == length(unique_dims),'unmatch number of cols');
            
            for i_dim = 1: length(unique_dims)
                
                idx_valid = hyper_pars.dim == unique_dims(i_dim);
                
                switch plotpars.err_metric
                    
                    case 'std'
                        
                        err_to_plot_test = plotpars.err_scale * err_test_cond_std;
                        err_to_plot_train = plotpars.err_scale * err_train_cond_std;
                        
                    case 'sem'
                        
                        err_to_plot_test = plotpars.err_scale * (1/sqrt(n_folds)).*err_test_cond_std;
                        err_to_plot_train = plotpars.err_scale * (1/sqrt(n_folds)).*err_train_cond_std;
                        
                end
                
                hh_test = ploterr(unique_lags, err_test_cond_mean(idx_valid,idx_cv_plot),[],err_to_plot_test(idx_valid,idx_cv_plot),'-','abshhy', 1.0);
                set(hh_test(1),'Color',plotpars.col_test(i_dim,:),'linewidth',plotpars.linewidth);
                set(hh_test(2),'Color',plotpars.col_test(i_dim,:),'linewidth',plotpars.bar_width);
                set(hh_test,'Parent',ah);
                
                hm_test = plot(ah,unique_lags, err_test_cond_mean(idx_valid,idx_cv_plot),'o');
                set(hm_test,'markeredgecolor',plotpars.col_test(i_dim,:) ,'markersize',plotpars.markersize,'markerfacecolor',[1 1 1]);
                
                if(plotpars.do_plot_train)
                    
                    hh_train = ploterr(unique_lags, err_train_cond_mean(idx_valid,idx_cv_plot),[],err_to_plot_train(idx_valid,idx_cv_plot),'-','abshhy', plotpars.bar_length);
                    set(hh_train(1),'Color',plotpars.col_train(i_dim,:),'linewidth',plotpars.linewidth);
                    set(hh_train(2),'Color',plotpars.col_train(i_dim,:),'linewidth',plotpars.bar_width);
                    
                    set(hh_train,'Parent',ah);
                    
                    hm_train = plot(ah,unique_lags, err_train_cond_mean(idx_valid,idx_cv_plot),'o');
                    set(hm_train,'markeredgecolor',plotpars.col_train(i_dim,:) ,'markersize',plotpars.markersize,'markerfacecolor',[1 1 1]);
                    
                end
                
                
            end
            
            
        case 'dim'
            
            assert(size(plotpars.col_test,1) == length(unique_lags),'unmatch number of cols');
            
            for i_dim = 1: length(unique_lags)
                
                idx_valid = hyper_pars.lag == unique_lags(i_dim);
                
                switch plotpars.err_metric
                    
                    case 'std'
                        
                        err_to_plot_test = plotpars.err_scale * err_test_cond_std;
                        err_to_plot_train = plotpars.err_scale * err_train_cond_std;
                        
                    case 'sem'
                        
                        err_to_plot_test = plotpars.err_scale * (1/sqrt(n_folds)).*err_test_cond_std;
                        err_to_plot_train = plotpars.err_scale * (1/sqrt(n_folds)).*err_train_cond_std;
                        
                end
                
                hh_test = ploterr(unique_dims, err_test_cond_mean(idx_valid,idx_cv_plot),[],err_to_plot_test(idx_valid,idx_cv_plot),'-','abshhy', plotpars.bar_length);
                set(hh_test(1),'Color',plotpars.col_test(i_dim,:),'linewidth',plotpars.linewidth);
                set(hh_test(2),'Color',plotpars.col_test(i_dim,:),'linewidth',plotpars.bar_width);
                set(hh_test,'Parent',ah);
                
                hm_test = plot(ah,unique_dims, err_test_cond_mean(idx_valid,idx_cv_plot),'o');
                set(hm_test,'markeredgecolor',plotpars.col_test(i_dim,:) ,'markersize',plotpars.markersize,'markerfacecolor',[1 1 1]);
                
                if(plotpars.do_plot_train)
                    
                    hh_train = ploterr(unique_dims, err_train_cond_mean(idx_valid,idx_cv_plot),[],err_to_plot_train(idx_valid,idx_cv_plot),'-','abshhy', plotpars.bar_length);
                    set(hh_train(1),'Color',plotpars.col_train(i_dim,:),'linewidth',plotpars.linewidth);
                    set(hh_train(2),'Color',plotpars.col_train(i_dim,:),'linewidth',plotpars.bar_width);
                    
                    set(hh_train,'Parent',ah);
                    
                    hm_train = plot(ah,unique_dims, err_train_cond_mean(idx_valid,idx_cv_plot),'o');
                    set(hm_train,'markeredgecolor',plotpars.col_train(i_dim,:) ,'markersize',plotpars.markersize,'markerfacecolor',[1 1 1]);
                    
                end
                
                
            end
   
    end
    
    xlabel(ah,plotpars.x_type);ylabel('cv-error');

    
end

h = gcf;
varargout{1} = h;
end