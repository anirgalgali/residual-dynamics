classdef ToyDynamicsModels
%{ This class is used to define the objects that correspond to the models 
%  of decisions and movements. See Figure 1.
% 
% Author : Aniruddh Galgali (Sep 2018)  
% 
%}
    
    properties
 
        model_type
        mu_init
        sigma_init
        noise_var
        sim_time
        ode_handle
        num_trajectories
        input_direction
        input_scale_factor
    end
    
properties(Constant)
        RAND_SEED = 541782;
        TIME_STEP = 0.001 ;
end

 
    
    methods
        
        function thisModel = ToyDynamicsModels(model,mu,sigma,noise_var,sim_time,...
                num_trajectories,ode_handle,input_direction,input_scale)
            % Constructor for the model object
            if(nargin == 9)
                thisModel.model_type = model;
                thisModel.mu_init = mu;
                thisModel.sigma_init= sigma;
                thisModel.noise_var = noise_var;
                thisModel.sim_time = sim_time;
                thisModel.num_trajectories = num_trajectories;
                thisModel.ode_handle = ode_handle;
                thisModel.input_direction = input_direction;
                thisModel.input_scale_factor = input_scale;
            end
        end
        
        function time_axes = generateTimeAxes(model)
           % Method to generate the time axis
           time_axes = 0: model.TIME_STEP: (model.sim_time - 1)*model.TIME_STEP;
        end
        
        function [x0] = generateInitialConditions(model,ic_type,varargin)
            % Method to generate initial conditions
            rng(model.RAND_SEED)
            switch ic_type
                
                case 'fixed'
                    assert(nargin > 2, 'not enough input arguments')                 
                    x0 = varargin{1};
                case 'random'
                    x0 = model.mu_init + chol(model.sigma_init)'*randn(2,model.num_trajectories);
            end
            
            
        end
        
        
        function [input] = generateInput(model,input_type,varargin)
             % Method to generate input drives. Input drives can either be
             % 'explicit', where the time-course of the input is
             % pre-specified, OR, 'optimized' where input is estimated to
             % so as the match the condition-average traajectories of a
             % reference model (Saddle point for models fo decisions and
             % Rotations for model of movement
            switch input_type
                
                case 'optimized'
                     
                    avgTraj = varargin{1};
                    nAvgs = size(avgTraj,3);
                    input_sign = varargin{2}; 
                    
                     for ii = 1:nAvgs
                         for tt = 1:size(avgTraj,2) - 1
                             
                             [x_dot_instantaneous] = model.ode_handle(avgTraj(:,tt,ii)');
                             input_avg(:,tt,ii) = (1./(model.TIME_STEP)).*(avgTraj(:,tt+1,ii) - avgTraj(:,tt,ii) - x_dot_instantaneous.*model.TIME_STEP) ;
     
                         end
                     end
    
                     input = zeros(2,model.sim_time - 1,model.num_trajectories);
                     input(:,:,input_sign == -1) = repmat(model.input_scale_factor.*input_avg(:,:,1),[1 1 sum(input_sign == -1)]);
                     input(:,:,input_sign == +1) = repmat(model.input_scale_factor.*input_avg(:,:,2),[1 1 sum(input_sign == +1)]);
                    
                case 'explicit-cnst'
                   
                    % explicitly specified input which is time-invariant
                    input_sign = varargin{1};
                    input = permute(repmat(model.input_scale_factor.*input_sign.*model.input_direction,[1 1 model.sim_time]),[1 3 2]);
                    
                   
                case 'explicit-timevarying'
                    % explicitly specified input which is time-varying
                    input_sign = varargin{1};
                    for tt = 1:model.sim_time
                        
                       input(:,tt,:) = model.input_scale_factor(:,tt).*input_sign.*model.input_direction;
                        
                    end
                    
                
            end
        end
        
        function [states] = generateStateTrajectories(model,x0,input,do_noisy,do_polar)
        % MEthod that genrates the latent state trajectories. Trajectories 
        % can be generated either using cartesian coordinates OR, using polar
        % coordinates (used only for models of movement)
            
            rng(model.RAND_SEED)

            for tt = 1:model.sim_time
                
                if(tt == 1)
                    

                    states(:,tt,:) = x0;
                    

                end
                
                if(do_polar)
                        
                       states_p(1,tt,:) = abs(squeeze(states(1,tt,:)) + 1i*squeeze(states(2,tt,:))); 
                       states_p(2,tt,:) = angle(squeeze(states(1,tt,:)) + 1i*squeeze(states(2,tt,:)));
                end
                
                
                if(tt < model.sim_time)
                    
                   
                   
                    
                    if(do_polar)
                        [s_dot_instantaneous] = model.ode_handle(squeeze(states_p(:,tt,:))');
                        
                        if(size(s_dot_instantaneous,1) == model.num_trajectories + 1 & s_dot_instantaneous(1) == 0)
                            
                           s_temp(1,:) = zeros(1,model.num_trajectories);
                           s_temp(2,:) = s_dot_instantaneous(2:end);
                           
                           s_dot_instantaneous = s_temp;
                        else
                           
                           s_dot_instantaneous = reshape(s_dot_instantaneous,[model.num_trajectories 2])';

                            
                            
                        end
                        
                        states_p(:,tt+1,:) = squeeze(states_p(:,tt,:)) + model.TIME_STEP.*(s_dot_instantaneous + squeeze(input(:,tt,:)));
                        states(1,tt+1,:) = squeeze(states_p(1,tt+1,:)).*cos(squeeze(states_p(2,tt+1,:)));
                        states(2,tt+1,:) = squeeze(states_p(1,tt+1,:)).*sin(squeeze(states_p(2,tt+1,:)));
                        
                    else
                        
                        [x_dot_instantaneous] = model.ode_handle(squeeze(states(:,tt,:))');
                        x_dot_instantaneous = reshape(x_dot_instantaneous,[model.num_trajectories 2])';
                        states(:,tt+1,:) = squeeze(states(:,tt,:)) + model.TIME_STEP.*(x_dot_instantaneous + squeeze(input(:,tt,:)));
                        
                        
                        
                    end
                    
                    
                    if(do_noisy)
                        states(:,tt+1,:) = squeeze(states(:,tt+1,:)) + chol(model.noise_var)'*randn(2,model.num_trajectories);
                    end
                    
                    
                end
                
            end


        end
        
        
        function [trajectories,flow] = generateGridEvolution(model,grid_points_X,grid_points_Y,time_steps,do_polar,forcing_term)
        % Method that unrolls the model trajectories for a few time-steps
        % starting from different state-space locations on a 2D grid (for
        % flow-field visualization purposes).
            
            
            [X1,X2] = meshgrid(grid_points_X,grid_points_Y);
            X_all = permute(cat(3,X1,X2),[3 1 2]);
            X_all = X_all(:,:);
            
            for ii = 1:size(X_all,2)
                
                for tt = 1:time_steps
                    
                    if(tt == 1)
                        
                        trajectories(1,tt,ii) = X_all(1,ii);
                        trajectories(2,tt,ii) = X_all(2,ii);
                        
                    end
                    
                    
                    if(do_polar)
                        
                        states_p(1,tt,ii) = abs(squeeze(trajectories(1,tt,ii)) + 1i*squeeze(trajectories(2,tt,ii)));
                        states_p(2,tt,ii) = angle(squeeze(trajectories(1,tt,ii)) + 1i*squeeze(trajectories(2,tt,ii)));
                        [s_dot_instantaneous] = model.ode_handle(squeeze(states_p(:,tt,ii))');
                         
                        if(tt < time_steps)
                            states_p(:,tt+1,ii) = squeeze(states_p(:,tt,ii)) + model.TIME_STEP.*(s_dot_instantaneous + forcing_term);
                            trajectories(1,tt+1,ii) = squeeze(states_p(1,tt+1,ii)).*cos(squeeze(states_p(2,tt+1,ii)));
                            trajectories(2,tt+1,ii) = squeeze(states_p(1,tt+1,ii)).*sin(squeeze(states_p(2,tt+1,ii)));
                        end
                        

                    
                    else
                        
                        if(tt < time_steps)
                            trajectories(:,tt+1,ii) = squeeze(trajectories(:,tt,ii)) + model.TIME_STEP.*(model.ode_handle(trajectories(:,tt,ii)'));
                        end
                    end
                end
                
                flow(:,ii) = (1./model.TIME_STEP).*(trajectories(:,2,ii) - trajectories(:,1,ii));
                
            end

        end
        
        
        function [norm_propagation,varargout] = compareTrueandEstimatedDynamics(model,A,initial_states,dirs_perturb,dirs_scale,time_steps,inputs)
            
            num_initial_states = size(initial_states,2);
            num_dirs_perturb = size(dirs_perturb,2);
            
            norm_ratio_true = zeros(num_dirs_perturb,num_initial_states);
            norm_ratio_pred = zeros(num_dirs_perturb,num_initial_states);
            res_true_all = zeros(2,time_steps,num_dirs_perturb,num_initial_states);
            res_pred_all = zeros(2,time_steps,num_dirs_perturb,num_initial_states);
            
            for ii = 1: num_initial_states
                
                x0 = initial_states(:,ii);
                idx_start = ii;
                
                for jj = 1:num_dirs_perturb
                    
                    x_const_true = zeros(2,time_steps);
                    x_perturb_true = zeros(2,time_steps);
                    x_perturb_pred = zeros(2,time_steps);
                    
                    for tt = 1:time_steps
                        
                        if(tt == 1)
                            
                            x_const_true(:,tt) = x0;
                            x_perturb_true(:,tt) = x0 + dirs_scale.*dirs_perturb(:,jj);
                            x_perturb_pred(:,tt) = x_perturb_true(:,tt) - x_const_true(:,tt);
                            
                        end
                        
                        if(tt < time_steps)
                            
                            x_dot_const = model.ode_handle(x_const_true(:,tt)');
                            x_dot_perturb = model.ode_handle(x_perturb_true(:,tt)');
                            
                            x_const_true(:,tt+1) = x_const_true(:,tt) + model.TIME_STEP.*(x_dot_const + inputs(:,idx_start + tt - 1));
                            x_perturb_true(:,tt+1) = x_perturb_true(:,tt) + model.TIME_STEP.*(x_dot_perturb + + inputs(:,idx_start + tt - 1));
                            
                            x_perturb_pred(:,tt+1) = A(:,:,idx_start + tt - 1)*x_perturb_pred(:,tt);
                            
                        end
                    end
                    
                    res_true = x_perturb_true - x_const_true;
                    
                    norm_ratio_true(jj,ii) = norm(res_true(:,end))./norm(res_true(:,1));
                    norm_ratio_pred(jj,ii) = norm(x_perturb_pred(:,end))./norm(x_perturb_pred(:,1));
                    
                    res_true_all(:,:,jj,ii) = res_true;
                    res_pred_all(:,:,jj,ii) = x_perturb_pred;
                    
                end
                
            end
            
            norm_propagation.true = norm_ratio_true;
            norm_propagation.pred = norm_ratio_pred;
            
            varargout{1} = res_true_all;
            varargout{2} = res_pred_all;
            
        end
        
    end
    
  
end