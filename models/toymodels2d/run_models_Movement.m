%{
This script generates the necessary data for the models of movement. We
simulate three different types of 2D dynamical models that produce average
trajectories with rotational structure
% 1) Autonomous rotational dynamics - in which dynamics are governed by a
     rotational dynamical system (without inputs)
% 2) Dynamic attractor (channel model) - in which the dynamics decay
     radially, but tend to push states through a "channel".
% 3) Point attractor - moving point attractor udner the influence of phasic
     inputs.
% Only model 3 has an input, whose time-course is inferred so as to obtain
a match with the condition-average trajectories of the rotations model
%
%
% Author: Aniruddh Galgali (Oct 2018)
%
%}
clearvars -except DIRS
clc

%% Initializing all the models

model_name = 'rotation';
rot_mat = @(theta) [cos(theta) sin(theta);-sin(theta) cos(theta)]  ;
pars_rotation.alpha = 1.2;
traj_rot_angle = -pi/120;
pars_rotation.A = [0 0;pars_rotation.alpha 0];
syms x y
x_dot = pars_rotation.A*[x;y];
model_ode = matlabFunction(x_dot, 'Vars', {[x, y]});
model1 = ToyDynamicsModels(model_name,[1;0],diag([0.01;1e-4]),diag([0.00003;0.00003]),...
    1000,4000,model_ode,[-1;1],[0;0]);    

model_name = 'channel';
pars_channel.alpha = 1.2;
pars_channel.beta = 20;
pars_channel.A = [-pars_channel.beta 0; pars_channel.alpha 0];
syms x y
x_dot = pars_channel.A*[x;y];
model_ode = matlabFunction(x_dot, 'Vars', {[x, y]});
model2 = ToyDynamicsModels(model_name,[1;0],diag([0.01;1e-4]),diag([0.0001;0.0001]),...
    1000,4000,model_ode,[pars_channel.beta;0],[1;1]);  
% 
model_name = 'moving_point_attractor';
pars_pointAtt.A = diag([-20; -40]);

syms x y
x_dot = pars_pointAtt.A*[x;y];
model_ode = matlabFunction(x_dot, 'Vars', {[x, y]});
model3 = ToyDynamicsModels(model_name,[1;0],diag([0.005;0.005]),diag([0.0003;0.0003]),...
    1000,4000,model_ode,[-1;1],[1;1]);

simPars.simulation_time  = generateTimeAxes(model1);

simPars.model1.pars = pars_rotation;
simPars.model1.model_obj = model1;

simPars.model2.pars = pars_channel;
simPars.model2.model_obj = model2;

simPars.model3.pars = pars_pointAtt;
simPars.model3.model_obj = model3;


%% Generate initial conditions for model

for ii = 1: model1.num_trajectories
    if(rand(1) <= 0.5)
        choiceLabels(ii) = -1;
    else
        choiceLabels(ii) = +1;
    end
    
end
mu_x0 = zeros(2,model1.num_trajectories);
mu_x0(:,choiceLabels == -1) = -1.*model1.mu_init.*ones(2,sum(choiceLabels == -1));
mu_x0(:,choiceLabels == +1) = model1.mu_init.*ones(2,sum(choiceLabels == +1));

for ii = 1:model1.num_trajectories
    
    x0_model1(:,ii) = mu_x0(:,ii) + chol(model1.sigma_init)'*randn(2,1);
end
%% Optimize/infer inputs so as to match condition-averages

[x0_model1] = generateInitialConditions(model1,'fixed',x0_model1);
[inputs_model1] = generateInput(model1,'explicit-cnst',choiceLabels);
[states_model1_noiseless] = generateStateTrajectories(model1,x0_model1,inputs_model1,false,true);
mean_model1(:,:,1) = mean(states_model1_noiseless(:,:,choiceLabels == -1),3);
mean_model1(:,:,2) = mean(states_model1_noiseless(:,:,choiceLabels == +1),3);
figure;plot(squeeze(mean_model1(1,:,1)),squeeze(mean_model1(2,:,1)),'-b')
hold on
plot(squeeze(mean_model1(1,:,2)),squeeze(mean_model1(2,:,2)),'-r')

[x0_model2] = generateInitialConditions(model2,'fixed',x0_model1);
[inputs_model2] = generateInput(model2,'explicit-cnst',ones(size(choiceLabels)));
[states_model2_noiseless] = generateStateTrajectories(model2,x0_model2,inputs_model2,false,true);
mean_model2(:,:,1) = mean(states_model2_noiseless(:,:,choiceLabels == -1),3);
mean_model2(:,:,2) = mean(states_model2_noiseless(:,:,choiceLabels == +1),3);
plot(squeeze(mean_model2(1,:,1)),squeeze(mean_model2(2,:,1)),'-.b')
hold on
plot(squeeze(mean_model2(1,:,2)),squeeze(mean_model2(2,:,2)),'-.r')


[x0_model3] = generateInitialConditions(model3,'fixed',x0_model1);
[inputs_model3] = generateInput(model3,'optimized',mean_model1,choiceLabels);
[states_model3_noiseless] = generateStateTrajectories(model3,x0_model3,inputs_model3,false,false);
mean_model3(:,:,1) = mean(states_model3_noiseless(:,:,choiceLabels == -1),3);
mean_model3(:,:,2) = mean(states_model3_noiseless(:,:,choiceLabels == +1),3);
plot(squeeze(mean_model3(1,:,1)),squeeze(mean_model3(2,:,1)),'-.b')
hold on
plot(squeeze(mean_model3(1,:,2)),squeeze(mean_model3(2,:,2)),'-.r')
%% Computing mean input for the models
unique_condition_labels = unique(choiceLabels);
for cidx = 1: length(unique_condition_labels)
    inputs.model1.mean(:,:,cidx) = mean(inputs_model1(:,:,choiceLabels == unique_condition_labels(cidx)),3);
    inputs.model1.choice_labels(cidx) = unique_condition_labels(cidx); 
    inputs.model1.time = simPars.simulation_time;
    inputs.model2.mean(:,:,cidx) = mean(inputs_model2(:,:,choiceLabels == unique_condition_labels(cidx)),3);
    inputs.model2.choice_labels(cidx) = unique_condition_labels(cidx);
    inputs.model2.time = simPars.simulation_time;
    inputs.model3.mean(:,:,cidx) = mean(inputs_model3(:,:,choiceLabels == unique_condition_labels(cidx)),3);
    inputs.model3.choice_labels(cidx) = unique_condition_labels(cidx);
    inputs.model3.time = simPars.simulation_time(1:end-1);
end
%% Generating noisy state trajectories an computing residuals

[states.model1.response] = generateStateTrajectories(model1,x0_model1,inputs_model1,true,true);
[states.model2.response] = generateStateTrajectories(model2,x0_model2,inputs_model2,true,true);
[states.model3.response] = generateStateTrajectories(model3,x0_model3,inputs_model3,true,false);
model_names = fieldnames(states);
unique_condition_labels = unique(choiceLabels);
for ii = 1: length(model_names)
   
    states.(model_names{ii}).choice_labels = choiceLabels;
    states.(model_names{ii}).time = simPars.simulation_time;
   
   for cidx = 1:length(unique_condition_labels)
      
       states.(model_names{ii}).condition_average(:,:,cidx) = mean(states.(model_names{ii}).response(:,:,choiceLabels ...
          == unique_condition_labels(cidx)),3);
      
       states.(model_names{ii}).residuals{cidx} = states.(model_names{ii}).response(:,:,choiceLabels == unique_condition_labels(cidx)) - ...
           repmat(states.(model_names{ii}).condition_average(:,:,cidx),[1 1 sum(choiceLabels == unique_condition_labels(cidx))]);
      
       [p,T,K] = size(states.(model_names{ii}).residuals{cidx});
       
       % adding a small amount of observation noise thats uncorrelated
       % across time makes the 2 stage regression more stable.
       obs_noise_small = reshape(chol(diag(1e-6.*ones(p,1)))'*randn(p,T*K),[p T K]);
       states.(model_names{ii}).residuals{cidx} = states.(model_names{ii}).residuals{cidx} + obs_noise_small;
   end
   
end
%% %% Fitting the model (No Cross-validation)
opt_alpha = 100;
opt_lag = 5;

fitPars.pastLag = opt_lag;
fitPars.regConst = opt_alpha;

for ii = 1: length(model_names)
    for cidx = 1:length(unique_condition_labels)
        
        time_labels = ones(length(states.(model_names{ii}).time),1);
        [analysis.(model_names{ii}).dynamics{cidx}.A,time_idxs_dyn] = estimateDynamicsTwoStageLS(states.(model_names{ii}).residuals{cidx}, ....
            opt_lag, opt_alpha, time_labels);
        
        analysis.(model_names{ii}).dynamics{cidx}.time = states.(model_names{ii}).time(time_idxs_dyn);
        
        [outA] = computeDerivedDynamicsQuantities(analysis.(model_names{ii}).dynamics{cidx}.A,'abs','right',time_labels(time_idxs_dyn),(states.(model_names{ii}).time(2) - states.(model_names{ii}).time(1)));
        outAnames = fieldnames(outA);
        
        for ff = 1: length(outAnames)
            
            analysis.(model_names{ii}).dynamics{cidx}.(outAnames{ff}) = outA.(outAnames{ff});
            
        end
        
        
        analysis.(model_names{ii}).dynamics{cidx}.fitPars = fitPars;
    end
end


%% Storing flow-field

grid_points = -1.5:0.1:1.5;
time_steps = 200;
for ii = 1: length(model_names)
    analysis_flow.(model_names{ii}).grid_points = grid_points;
    analysis_flow.(model_names{ii}).time_steps = time_steps;
    if(ii == 1)
        [analysis_flow.(model_names{ii}).trajectories,analysis_flow.(model_names{ii}).flow] = generateGridEvolution(simPars.(model_names{ii}).model_obj,grid_points,time_steps,true,[0;0]);
    elseif(ii == 3)
        [analysis_flow.(model_names{ii}).trajectories,analysis_flow.(model_names{ii}).flow] = generateGridEvolution(simPars.(model_names{ii}).model_obj,grid_points,time_steps,false,[0;0]);
    else
        [analysis_flow.(model_names{ii}).trajectories,analysis_flow.(model_names{ii}).flow] = generateGridEvolution(simPars.(model_names{ii}).model_obj,grid_points,time_steps,true,[pars_channel.beta;0]);
    end
end

%%

movement_models.simPars = simPars;
movement_models.analysis = analysis;
movement_models.data.states = states;
movement_models.data.inputs = inputs;
movement_models.data.phase = analysis_flow;

if(~isunix)
    save_path = '';
else
    save_path = fullfile(DIRS.analysis,'/simulations/toymodels/');
end
save(fullfile(save_path,strcat('movement_models_newAnalysisPipeline_',date,'.mat')),'movement_models');
