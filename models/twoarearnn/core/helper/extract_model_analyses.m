function [result, self_rec, fwd, fbk, varargout] = extract_model_analyses(file_path,file_name,model_pars)
%{ A helper function that extracts the results of the two-area network 
%  simulations stored to disk

% INPUT
% - file_path (string) - path to stored file
% - file_name (string) - name of .mat file with stored models 
% - model_pars (table) - of size (n_models x 3) - A table specifying a specifc
%   model that you may want to extract for further analysis. Each row of the 
%   table specifies a model (in total n_models model), and columns indicate
%   the parameters (self_recurrence, feed-fwd and feedback stengths)
%
% OUTPUT
% - result (struct) - The full data structure containing the result for all
% simulated models
% - self_rec (array - n_cross x n_self) - a matrix specifing the self recurrence strengths for
% all simulated models, where n_cross is the total number of swept fwd/fbk
% weights
% - fwd (array - n_cross x n_self) - a matrix specifying the feed-fwd strengths of all
% simulated models
% - fbk (array) - a matrix specifying the feedback strengths of all
% simulated models
% - varargout - (i) model_out (cell_array) - of size n_models x 1 containing 
% the results for each specified model in model_pars.
%
% Auhor: Aniruddh Galgali
%}

%% Loading stored results

% Load file from disk
R = load(fullfile(file_path,file_name)); % loads a variable called 'result'
result = R.result;

% Extracting the network parameters for all models
self_rec = cell2mat(arrayfun(@(x) x.self_ppc,result.net_pars,'uni',false));
fwd = cell2mat(arrayfun(@(x) x.fwd_ppc_pfc,  result.net_pars,'uni',false));
fbk = cell2mat(arrayfun(@(x) x.fbk_pfc_ppc,  result.net_pars,'uni',false));
c_count = 1;
out = {};
% Extracting only those models specified in model_pars
for imdl = 1: size(model_pars,1)
    [r,c] = find((abs(self_rec - model_pars.self(imdl)) < 1e-3) & ...
        (abs(fwd - model_pars.fwd(imdl)) < 1e-3) & ...
        (abs(fbk - model_pars.fbk(imdl)) < 1e-3));
    
    if(isempty(r) || isempty(c))
        
        continue;
        
    else
        
        out{c_count}.result_cv = squeeze(result.analysis(r,c,:));
        out{c_count}.data = squeeze(result.data{r,c});
        out{c_count}.net_pars = result.net_pars(r,c);
        out{c_count}.input_pars = result.input_pars;
        out{c_count}.obs_pars = result.obs_pars;
        out{c_count}.sim_pars = result.sim_pars;
        out{c_count}.areas = result.areas;
        c_count = c_count + 1;
        
    end
    
end

varargout{1} = out;

end