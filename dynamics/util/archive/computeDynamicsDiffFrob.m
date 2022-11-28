function [frobDiff] = computeDynamicsDiffFrob(A1,varargin)

if(nargin == 1)
    A2 = A1;
elseif(nargin == 2)
    A2 = varargin{1};
end

nTime = size(A1,3);
frobDiff = NaN(nTime,nTime); 

for tt1 = 1:nTime
    for tt2 = 1:nTime
        frobDiff(tt1,tt2) = norm(A1(:,:,tt1) - A2(:,:,tt2),'fro');
    end
    
end
end