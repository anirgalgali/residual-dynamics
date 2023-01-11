%{
This script generates the necessary data for the models of decisions. We
simulate three different types of 2D dynamical models that integrate an
input towards choice
% 1) Saddle point model (unstable integration)
% 2) Line attractor model (perfect integration)
% 3) Point attractor (leaky integration)
% All 3 models are driven by external inoputs, but only the inputs for the
saddle model are specified explicityl, whereas the inputs for the ther two
models are inferred so as to have a match in the condition-averages across
all 3 model types
%
%
% Author: Aniruddh Galgali (Oct 2018)
%
%}
clearvars -excet DIRS
clc

%% Initializing all the models

model_name = 'simple_saddle';
pars_saddle_ode.a = -6;
pars_saddle_ode.b = 20;
pars_saddle_ode.c = -10;
[model_saddle] = generateSaddleDynamics(pars_saddle_ode,model_name);
inputs_to_saddle = [repmat([8;-8],[1 1000])];
model1 = ToyDynamicsModels(model_name,[0;0],diag([0.005;0.005]),diag([0.0001;0.0001]),...
    1000,4000,model_saddle.model_ode,[-1;1],inputs_to_saddle);   

model_name = 'non_normal';
pars_nonNormal_la.l0 = [-1;1];
pars_nonNormal_la.r0 = [1;0];
pars_nonNormal_la.dyn_eigenval = [0 -5];
pars_nonNormal_la.A = generateParametersforNonNormalDynamics(pars_nonNormal_la);

syms x y
x_dot = pars_nonNormal_la.A*[x;y];
model_ode = matlabFunction(x_dot, 'Vars', {[x, y]});
model2 = ToyDynamicsModels(model_name,[0;0],diag([0.005;0.005]),diag([0.00003;0.00003]),...
    1000,4000,model_ode,[-1;1],[1;1]);  


model_name = 'moving_point_attractor';
pars_pointAtt.A = diag([-20; -40]);

syms x y
x_dot = pars_pointAtt.A*[x;y];
model_ode = matlabFunction(x_dot, 'Vars', {[x, y]});
model3 = ToyDynamicsModels(model_name,[0;0],diag([0.005;0.005]),diag([0.0003;0.0003]),...
    1000,4000,model_ode,[-1;1],[1;1]);

simPars.simulation_time  = generateTimeAxes(model1);

simPars.model1.pars = pars_saddle_ode;
simPars.model1.fixed_points = model_saddle.fixed_points;
simPars.model1.model_jacobian = model_saddle.model_jac;
simPars.model1.model_obj = model1;

simPars.model2.pars = pars_nonNormal_la;
simPars.model2.model_obj = model2;

simPars.model3.pars = pars_pointAtt;
simPars.model3.model_obj = model3;

%% Generate initial conditions, optimize/infer inputs so as to match condition-averages

[x0_model1] = generateInitialConditions(model1,'random');
choiceLabels = zeros(1,model1.num_trajectories);
for ii = 1: model1.num_trajectories
    if(rand(1) <= 0.5)
        choiceLabels(ii) = -1;
    else
        choiceLabels(ii) = +1;
    end
    
end
[inputs_model1] = generateInput(model1,'explicit-timevarying',choiceLabels);
[states_model1_noiseless] = generateStateTrajectories(model1,x0_model1,inputs_model1,false,false);
mean_model1(:,:,1) = mean(states_model1_noiseless(:,:,choiceLabels == -1),3);
mean_model1(:,:,2) = mean(states_model1_noiseless(:,:,choiceLabels == +1),3);

[x0_model2] = generateInitialConditions(model2,'fixed',x0_model1);
[inputs_model2] = generateInput(model2,'optimized',mean_model1,choiceLabels);

[x0_model3] = generateInitialConditions(model3,'fixed',x0_model1);
[inputs_model3] = generateInput(model3,'optimized',mean_model1,choiceLabels);

%% Computing mean input for the models

unique_condition_labels = unique(choiceLabels);
for cidx = 1: length(unique_condition_labels)
    inputs.model1.mean(:,:,cidx) = mean(inputs_model1(:,:,choiceLabels == unique_condition_labels(cidx)),3);
    inputs.model1.choice_labels(cidx) = unique_condition_labels(cidx); 
    inputs.model1.time = simPars.simulation_time;
    inputs.model2.mean(:,:,cidx) = mean(inputs_model2(:,:,choiceLabels == unique_condition_labels(cidx)),3);
    inputs.model2.choice_labels(cidx) = unique_condition_labels(cidx);
    inputs.model2.time = simPars.simulation_time(1:end-1);
    inputs.model3.mean(:,:,cidx) = mean(inputs_model3(:,:,choiceLabels == unique_condition_labels(cidx)),3);
    inputs.model3.choice_labels(cidx) = unique_condition_labels(cidx);
    inputs.model3.time = simPars.simulation_time(1:end-1);
end

%% Generating noisy state trajectories an computing residuals

[states.model1.response] = generateStateTrajectories(model1,x0_model1,inputs_model1,true,false);
[states.model2.response] = generateStateTrajectories(model2,x0_model2,inputs_model2,true,false);
[states.model3.response] = generateStateTrajectories(model3,x0_model3,inputs_model3,true,false);
unique_condition_labels = unique(choiceLabels);

model_names = fieldnames(states);

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

%% Fitting the model (No Cross-validation)

opt_alpha = 100;
opt_lag = 5;

fitPars.pastLag = opt_lag;
fitPars.regConst = opt_alpha;

for ii = 1: length(model_names)
    for cidx = 1:length(unique_condition_labels)
        
        time_labels = ones(length(states.(model_names{ii}).time),1);
        [analysis.(model_names{ii}).dynamics{cidx}.A,time_idxs_dyn] = ...
            estimateDynamicsTwoStageLS(states.(model_names{ii}).residuals{cidx}, ...
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
grid_points = -3:0.1:3;
time_steps = 300;
for ii = 1: length(model_names)
   analysis_flow.(model_names{ii}).grid_points = grid_points;
   analysis_flow.(model_names{ii}).time_steps = time_steps;
   [analysis_flow.(model_names{ii}).trajectories,analysis_flow.(model_names{ii}).flow] = generateGridEvolution(simPars.(model_names{ii}).model_obj,grid_points,grid_points,time_steps,false);
    
end

%% Saving all results to disk

integration_models.simPars = simPars;
integration_models.analysis = analysis;
integration_models.data.states = states;
integration_models.data.inputs = inputs;
integration_models.data.phase = analysis_flow;
if(~isunix)
    save_path = '';
else
    save_path = fullfile(DIRS.analysis,'/simulations/toymodels/');
end
save(fullfile(save_path,strcat('integration_models_newAnalysisPipeline_',date,'.mat')),'integration_models');

