function [H] = plotStateTrajectories(X,pars,labels,varargin)
%{ This function plots the trajectories in a 2D state-space.
% Inputs
% - X (3D tensor) - either of size n_dim x n_times x n_trials/n_conds. 
%   Specifies the data to be plotted 
% - pars (struct) - high-level figure properties
% - labels (array) - of size n_trials x 1 or n_conds x 1, which indicates 
%   the condition label associated with each trial/condition
% - varargin
%   - h - figure handle
%   - ah - axis handle
%}
if(nargin > 3)
   h = varargin{1};
   ah = varargin{2};
%    hold(ah);
elseif(nargin == 3)
    
    figure; ah = gca; hold(ah,'plotbox',[1 1 1]);
    h = gcf;
end
if(isfield(pars,'figPosition'))
    set(h,'Position',pars.figPosition);
end
% Selecting a random subset of trials
if(isfield(pars,'num_trials_to_plot'))
    rng(pars.rseed);
    trial_idxs = randperm(size(X,3),pars.num_trials_to_plot);
    X = X(:,:,trial_idxs);
    labels = labels(trial_idxs);
    
end

if(size(pars.Color,1) == 1)

    line_cols = repmat(pars.Color,[length(labels) 1]);

elseif(size(pars.Color,1) == 2)
   
    line_cols = NaN(length(labels),3);
    unique_labels = unique(labels);
    
    for cc = 1: length(unique_labels)
       
        line_cols(labels == unique_labels(cc),:) = ...
            repmat(pars.Color(cc,:),[sum(labels == unique_labels(cc)) 1]);  
        
    end
    
else
    
    error('excess conditions');
    
end



h_l = plot(ah,squeeze(X(1,:,:)),squeeze(X(2,:,:)),'linewidth',pars.lineWidth);
set(h_l, {'color'}, num2cell(line_cols,2))

if(pars.doMarkers)
%     
    X_m = X(:,1:pars.step:end,:);
    
    if(size(pars.markerEdgeCol,1) == 1)

        marker_edge_cols = repmat(pars.markerEdgeCol,[length(labels) 1]);

    elseif(size(pars.markerEdgeCol,1) == 2)
   
        marker_edge_cols = NaN(length(labels),3);
        unique_labels = unique(labels);

        for cc = 1: length(unique_labels)

            marker_edge_cols (labels == unique_labels(cc),:) = ...
                repmat(pars.markerEdgeCol(cc,:),[sum(labels == unique_labels(cc)) 1]);  

        end
    
    else
    
        error('excess conditions');
    
    end
    
    
    if(size(pars.markerFaceCol,1) == 1)

        marker_face_cols = repmat(pars.markerFaceCol,[length(labels) 1]);

    elseif(size(pars.markerFaceCol,1) == 2)
   
        marker_face_cols = NaN(length(labels),3);
        unique_labels = unique(labels);

        for cc = 1: length(unique_labels)

            marker_face_cols(labels == unique_labels(cc),:) = ...
                repmat(pars.markerFaceCol(cc,:),[sum(labels == unique_labels(cc)) 1]);  

        end
    
    else
    
        error('excess conditions');
    
    end
    
    for jj = 1:size(marker_edge_cols,1)
%         
        h_m = plot(ah,squeeze(X_m(1,:,jj)),squeeze(X_m(2,:,jj)),'o','markersize',pars.markersize);
        set(h_m, 'MarkerFaceColor', marker_face_cols(jj,:),'MarkerEdgeColor', marker_edge_cols(jj,:))
%         
    end
%    
%    
end

H.fig_handle = h;
H.axes_handle = ah;

end