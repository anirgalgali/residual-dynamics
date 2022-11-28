function [theta] = subspacea(F,G,A)
%SUBSPACEA angles between subspaces
%  subspacea(F,G,A)
%  Finds all min(size(orth(F),2),size(orth(G),2)) principal angles
%  between two subspaces spanned by the columns of matrices F and G 
%  in the A-based scalar product x'*A*y, where A
%  is Hermitian and positive definite. 
%  COS of principal angles is called canonical correlations in statistics.  
%  [theta,U,V] = subspacea(F,G,A) also computes left and right
%  principal (canonical) vectors - columns of U and V, respectively.
%
%  If F and G are vectors of unit length and A=I, 
%  the angle is ACOS(F'*G) in exact arithmetic. 
%  If A is not provided as a third argument, than A=I and 
%  the function gives the same largest angle as SUBSPACE.m by Andrew Knyazev,
%  see
%  http://www.mathworks.com/matlabcentral/fileexchange/Files.jsp?type=category&id=&fileId=54
%  MATLAB's SUBSPACE.m function is still badly designed and fails to compute 
%  some angles accurately.
%
%  The optional parameter A is a Hermitian and positive definite matrix,
%  or a corresponding function. When A is a function, it must accept a
%  matrix as an argument. 
%  This code requires ORTHA.m, Revision 1.5.8 or above,
%  which is included. The standard MATLAB version of ORTH.m
%  is used for orthonormalization, but could be replaced by QR.m.
%  
%  Examples: 
%  F=rand(10,4); G=randn(10,6); theta = subspacea(F,G);
%  computes 4 angles between F and G, while in addition 
%  A=hilb(10); [theta,U,V] = subspacea(F,G,A);
%  computes angles relative to A and corresponding vectors U and V. 
%  
%  The algorithm is described in A. V. Knyazev and M. E. Argentati,
%  Principal Angles between Subspaces in an A-Based Scalar Product: 
%  Algorithms and Perturbation Estimates. SIAM Journal on Scientific Computing, 
%  23 (2002), no. 6, 2009-2041.
%  http://epubs.siam.org/sam-bin/dbq/article/37733
%  Tested under MATLAB R10-14
%  Copyright (c) 2000 Andrew Knyazev, Rico Argentati
%  Contact email: knyazev@na-net.ornl.gov
%  License: free software (BSD)
%  $Revision: 4.5 $  $Date: 2005/6/27


% Edited by Aniruddh Galgali : 11/2020
% Edit includes a reordering of the values such that they are consistent with the directions specified in F. 
threshold=sqrt(2)/2; % Define threshold for determining when an angle is small
if size(F,1) ~= size(G,1)
   subspaceaError(['The row dimension ' int2str(size(F,1)) ...
         ' of the matrix F is not the same as ' int2str(size(G,1)) ...
         ' the row dimension of G'])
end
if nargin<3  % Compute angles using standard inner product
   
   % Trivial column scaling first, if ORTH.m is used later 
   for i=1:size(F,2),
     normi=norm(F(:,i),inf);
     %Adjustment makes tol consistent with experimental results
     if normi > eps^.981
       F(:,i)=F(:,i)/normi;
       % Else orth will take care of this
     end
   end
   for i=1:size(G,2),
     normi=norm(G(:,i),inf);
     %Adjustment makes tol consistent with experimental results
     if normi > eps^.981
       G(:,i)=G(:,i)/normi;
       % Else orth will take care of this
     end
   end
  % Compute angle using standard inner product
  
  QF = orth(F);      %This can also be done using QR.m, in which case
  QG = orth(G);      %the column scaling above is not needed 
  
  
  q = min(size(QF,2),size(QG,2));
  [Ys,s,Zs] = svd(QF'*QG,0);
  if size(s,1)==1
    % make sure s is column for output
    s=s(1);
  end
  s = min(diag(s),1);
  
  
  if(length(s) > 1)
      
      dp = diag(real(abs(QF'*F)));
      if(abs(sum(dp) - 0) <= 1e-5)
          s = flipud(s);
          
      end
      
  end
      
  theta = max(acos(s),0);
end