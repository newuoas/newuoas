function [ p, itn ] = mss( g, delta, tolmin, update_flag, S, Y, a, b, RHO, gamma_inv );
% 
%  Implementation of the MSS method to solve
%
%          min g'p + \half p'Bp   s.t. ||p|| <= delta
%
%  where B is an L-BFGS matrix with B_0 = (gamma_inv)*I 
%
%  More information on the MSS method can be found in:
%
%     "MSS: MATLAB Software for L-BFGS Trust-Region 
%       Subproblems for Large-Scale Optimization"
%
%  Authors: Jennifer Erway and Roummel Marcia
%
% The technical report and software are available at 
%   www.wfu.edu/~erwayjb/publications.html
%   www.wfu.edu/~erwayjb/software.html
%
%
%  Copyright (2013): Jennifer Erway and Roummel Marcia
%
% This code is distributed under the terms of the GNU General Public
% License
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
%---------------------------------------------------------------------------%
%
% Inputs: 
% =======
%   * g, delta -- as defined above
%   * tolmin -- tolerance for the More'-Sorensen method
%   * update_flag == 0 if the BFGS pairs (S,Y) have not been updated
%                    since last called; else, update_flag==1
%                    When used inside a trust-region method, this saves
%                    recomputing the variables CTC0 and C when the pairs
%                    (S,Y) have not been updated.  
%                    When used outside a trust-region algorithm setting,
%                    set update_flag=1.
%
%
% Global variables: 
% =================
%  * a,b are the rank 1 updates computed according to unrolling.m (see
%                website above or Procedure 7.6 in J. Nocedal,
%                S.J. Wright, Numerical Optimization, Springer-Verlag, New
%                York, second edition, 2006.)
%  * gamma_inv -- defined above
%  * m         -- the current number of limited memory updates
%  * (S,Y)     -- the bfgs pairs (most recent pair stored in S(:,m), Y(:,m))
%  * RHO(i)    --- holds 1/(S(:,i)'*Y(:,i))
%  * CTC0 and C are defined below (see "update_flag" above)

%global C CTC0;

itn     = 0;
n       = size( g, 1 );
m       = size( S, 2 );
max_itn = min(n,100);
sigma   = 0;
epsil   = sqrt(eps);

if (update_flag==1)
  C = zeros(n,2*m);
  C(:,[1:2:2*m]) = a;
  C(:,[2:2:2*m]) = b;
  CTC0 = C'*C; 
end

if m > 0
  gamma_qN = 1/gamma_inv;
  p        = two_loop(-g,gamma_qN,m, S, Y, RHO );
  pnorm    = norm(p);

%  Already false !!

  if ((pnorm<=delta) || (abs(pnorm-delta) <= tolmin*delta ))  
    %done
		
  else %step 4--6 
    CTg = C'*g;
    for itn=1:max_itn
      % update sigma
      phi  = 1/pnorm - 1/delta;
      if sigma<epsil
	gradp = two_loop(-p,gamma_qN,m, S, Y, RHO ); 
	sigma = 0;
      else
	gradp = newrec(gamma_inv,sigma,-p,CTC0,-C'*p,C);  
      end
      phi_prime = -(p'*gradp)/((pnorm)^3);
      sigma     =  sigma - phi/phi_prime;
      
      if sigma < epsil
	p     = two_loop(-g,gamma_qN,m, S, Y, RHO);
	sigma = 0;
      else
	p     = newrec(gamma_inv,sigma,-g,CTC0,-CTg,C);  
      end
      
      pnorm = norm(p);
      itn   = itn+1;  
      if (abs(pnorm-delta) <= tolmin*delta )
	break;  
      end
    end %for loop
  end %if loop
  
else

  p        = -g;
  pnorm    = norm(p);
  
  if ((pnorm<=delta) || (abs(pnorm-delta) <= tolmin*delta ))  
    %done
    
  else %step 4--6 
    for itn=1:max_itn
      % update sigma
      phi  = 1/pnorm - 1/delta;
      if sigma < epsil
	gradp  = -p;
	sigma  =  0;
      else 
	gradp  = -(1/(1+sigma))*p;
      end
      phi_prime = -(p'*gradp)/((pnorm)^3);
      sigma     =  sigma - phi/phi_prime;

      if sigma < epsil
	p     = -g;
	sigma =  0;
      else
	p = -(1/(1+sigma))*g;
      end
      
      pnorm = norm(p);
      itn   = itn+1;  
      if (abs(pnorm-delta) <= tolmin*delta )
	break;  
      end
    end %for loop
  end %if loop
  
end

if itn>=max_itn
  fprintf( 'mss could not find a solution\n' );
%  keyboard
end


