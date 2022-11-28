function [H_t,varargout] = computeTimeVaryingHankelMatrix(X, time, q, varargin)
%{ 

Computes a sequence of time-varying hankel matrices H(t) using the data in 
the 3-dimensional tensor X. It automatically splits the computation of the 
hankel matrices across different event alingments of the data according to
the event label associated with each time-point in X. These event labels
are specified as an optional argument. 

Inputs : X - 3D data tensor of size n_dim x n_time x n_trials
         time - time stamps associated with the data in X. An array of size
                1 x n_time
         q - specifies the hankel order i.e the number of lags of X to used
             to compute H(t). The number of "forward" and "backward" lags 
             are tied to one another (scalar)

         OPTIONAL (in-order):
            1) "time-labels" - A label associated with each time point 
                indicating the event alignment associated with that time
                bin. An array of size 1 x n_time

Outputs: H_t - A sequence of time-varying hankel matrices H(t), which is a 
               3D tensor of size (n_dim * q) x (n_dim * q) x n_times_hankel
               , where n_times_hankel = ((sum_{i} n_times(i)) - n_align*(2*q + 1)
               , where n_times(i) s the number of time-bins in the ith alignment 
               and n_align is the total number of event alignment, which is
               given by the number of unique elements of time-labels.

               OPTIONAL (in-order):
               1) time_all - An array of size n_times_hankel x 1, that 
                  specifies the time stamp associated with the corresponding 
                  H_t
               2) time_all_label - An array of size n_times_hankel x 1, that 
                  specifies the event alignment associated with each time-bin 
                  in time_all

        Author : Aniruddh Galgali, Dec 2017

%} 

if(nargin == 3)
    time_labels = ones(length(time),1);
elseif(nargin > 3)
     time_labels = varargin{1};
end
unique_time_labels = unique(time_labels);
n_align = length(unique_time_labels);

[n_dim, ~, ~] = size(X);

n_times_in_alignment = NaN(n_align , 1);
for ialign = 1:n_align
   idx_times = time_labels== unique_time_labels(ialign);
   n_times_in_alignment(ialign) = sum(idx_times);  
end

H_t = NaN(n_dim*q, n_dim*q, sum(n_times_in_alignment)- n_align*(2*q + 1));
time_all = [];
time_all_label = [];
t_count = 1;

for ialign = 1:n_align
    
    idx_times = time_labels== unique_time_labels(ialign);
    X_seg = X(:,idx_times,:);
    time_seg = time(idx_times);
    time_seg_label = time_labels(idx_times);
    X_seg = permute(X_seg,[3 1 2]);
        
    for tt = 1: length(time_seg) - 2*q + 1
       
        X_past = X_seg(:,:, tt + q -1 : - 1: tt);
        X_past = X_past(:,:);

        X_future = X_seg(:,:,tt + q : tt + 2*q  -1);
        X_future = X_future(:,:);

        H_t(:,:,t_count) = (1./(size(X_seg,1) - 1)).*(X_future'*X_past);
        
        t_count = t_count + 1;
        
    end

    time_seg = time_seg(q + 1: end - q + 1);
    time_seg_label = time_seg_label(q + 1: end - q + 1);
    time_all = [time_all time_seg];
    time_all_label = [time_all_label time_seg_label];
end

varargout{1} = time_all;
varargout{2} = time_all_label;

end