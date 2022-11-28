function pvar = percvar(data,pred)
% percvar percentage of variance explained
%
% 2003-09-10 VM made it
%
% pvar = percvar(data,pred)

if size(data) == size(pred)
   exp = data;
   mod = pred;
elseif size(data,1) > size(pred,1)
   exp = data;
   mod = repmat(pred,size(exp,1),1);
elseif size(data,2) > size(pred,2)
   exp = data;
   mod = repmat(pred,1,size(exp,2));
else
   error('bad input size');
end   

pvar = 100*(1 - sum((exp(:) - mod(:)).^2)/sum((exp(:) - mean(exp(:))).^2));

