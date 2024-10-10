function [out] = newrec(b0,sigma,z,CTC0,CTz0,C)
%
%
% This function directly solves systems of the form (B+\sigma I)x=z, 
% where B is an L-BFGS quasi-Newton matrix.  
%
% "Limited-Memory BFGS Systems with Diagonal Updates" 
% by Jennifer Erway and Roummel Marcia
%
% Copyright (2011): Jennifer Erway and Roummel Marcia
%
% The technical report and software are available at 
% www.wfu.edu/~erwayjb/publications.html
% www.wfu.edu/~erwayjb/software.html
%
% 
% This code is distributed under the terms of the GNU General Public License
% 2.0.
%
% Permission to use, copy, modify, and distribute this software for
% any purpose without fee is hereby granted, provided that this entire 
% notice is included in all copies of any software which is or includes
% a copy or modification of this software and in all copies of the
% supporting documentation for such software.
% This software is being provided "as is", without any express or
% implied warranty.  In particular, the authors do not make any
% representation or warranty of any kind concerning the merchantability
% of this software or its fitness for any particular purpose.
%---------------------------------------------------------------------------


% Inputs: a(:,k) and b(:,k) make up the kth L-bfgs update and
% B_0 = b0*I is the initial L-BFGS matrix, i.e., 
% B_k = B_0 - sum_{i=1}^k a(:,i)a(:,i)^T + sum_{i=1}^{k} b(:,i)b(:,i)^T

% initialization
k      = size(C,2);
inv_c0 = 1/( b0 + sigma );

PTC  = zeros(k,k);
PTz  = zeros(k,1);
v    = zeros(k,1);
p    = zeros(size(C,1),k);

CTC  = CTC0/(b0+sigma); 
CTz  = CTz0/(b0+sigma);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recursion formula
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j = 1:k;
    for i = 1:k
       PTC(j,i) = CTC(j,i);
       for t = 1:j-1
           PTC(j,i) = PTC(j,i) + ((-1)^(t+1))*v(t)*PTC(t,j)*PTC(t,i);
       end
    end

    PTz(j,1) = CTz(j,1);
    for t = 1:j-1
        PTz(j,1) = PTz(j,1) + ((-1)^(t+1))*v(t)*PTC(t,j)*PTz(t,1);
    end
    v(j) = 1/(1 + ((-1)^j)*PTC(j,j));
end

invCz  = inv_c0*z;
for j = 1:k;
    p(:,j) = inv_c0*C(:,j);
    for t = 1:j-1
        p(:,j) = p(:,j) + ((-1)^(t+1)*v(t)*(PTC(t,j)))*p(:,t);
    end
    
    % compute inner products using precomputed values in CTC
    invCz = invCz + (-1)^(j+1)*v(j)*PTz(j,1)*p(:,j);

end

out = invCz;



