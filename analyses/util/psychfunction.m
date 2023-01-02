function [xx,yy,nn] = psychfunction(xvalues,ichoice,itrial,varargin)
% PSYCHFUNCTION compute psychometric function
%
% Inputs:
% xvalues: abscissa values (1 x ntrials)
% ichoice: indexes of trials where the option was chosen (can be logical)
% itrial: subset of trials to consider (can be logical)
% varargin: parameters to compute itrial inside the function
%  consists of sequences of variables like
%  (...,'index',isub,...) or
%  (...,'value',subvaluelist,subvalue,...)
%
% Outputs:
% xx: absissa
% yy: ordinate (percentages)
% nn: ordinate (total number of choices)
%
% [xx,yy,nn] = psychfunction(xvalues,ichoice,itrial,varargin)

% The trials to consider
if nargin <3 || isempty(itrial)
   itrial = 1:length(xvalues);
end

% Convert to indeces if necessary
if islogical(ichoice)
    ichoice = find(ichoice);
end
if islogical(itrial)
    itrial = find(itrial);
end

% The trials to look at
ptrial = itrial;
ninput = length(varargin);
inext = 1;
while inext <= ninput
   switch varargin{inext}
      case 'index'
          if islogical(varargin{inext+1})
              ptrial = intersect(ptrial,find(varargin{inext+1}));
          else
              ptrial = intersect(ptrial,varargin{inext+1});
          end
         inext = inext+2;
      case 'value'
         ptrial = intersect(ptrial,find(varargin{inext+1} == varargin{inext+2}));
         inext = inext+3;
      otherwise
         error('Don''t understand these inputs');
   end
end

% All the x-axis values
xx = unique(xvalues(ptrial));

% Get the corresponding choices
yy = zeros(size(xx));
nn = zeros(size(xx));

for ix = 1:length(xx)
   ipass = intersect(ptrial,find(xvalues == xx(ix)));
   yy(ix) = length(intersect(ichoice,ipass))/length(ipass);
   nn(ix) = length(ipass);
end

return

