function [ x, error, iter, flag, prered ] = ...
         pcg( matvec, b, precond, max_it, tol, radius, reorth, S, Y )

%  -- Iterative template routine --
%     Univ. of Tennessee and Oak Ridge National Laboratory
%     October 1, 1993
%     Details of this algorithm are described in "Templates for the
%     Solution of Linear Systems: Building Blocks for Iterative
%     Methods", Barrett, Berry, Chan, Demmel, Donato, Dongarra,
%     Eijkhout, Pozo, Romine, and van der Vorst, SIAM Publications,
%     1993. (ftp netlib2.cs.utk.edu; cd linalg; get templates.ps).
%
%  [x, error, iter, flag] = cg( matvec, x, b, dim, precond, max_it, tol)
%
% cg.m solves the symmetric positive definite linear system Ax=b
% using the Conjugate Gradient method with preconditioning.
%
% input   matvec   the handle of a function computing the
%                  product of A times a vector.
%         b        REAL right hand side vector
%         precond  the handle of a function computing M\x
%         max_it   INTEGER maximum number of iterations
%         tol      REAL relative error tolerance
%
% output  x        REAL solution vector
%         error    REAL relative error norm
%         iter     INTEGER number of iterations performed
%         flag     INTEGER: 0 = solution found to tolerance
%                           1 = no convergence given max_it
%
% Updated August 2006; rbarrett@ornl.gov. (See ChangeLog for details.)
% Modified to cope with function handles for the product with A and M
% (Ph. Toint, 9 V 2016) and for a simple truncated CG trust-region
% mechanism (Ph. Toint, 25 VI 2017).
%
% =============================================================================

  TR_use  = 1;
  verbose = 0;
  
% ---------------
% Initializations
% ---------------

  if ( isnumeric( precond ) && size( precond, 1 ) == 1 && precond <= 0 )
     use_precond = 0;
  else
     use_precond = 1;
  end
  
  flag   = 0;
  iter   = 0;
  prered = 0;
  alpha  = 0.0;
  beta   = 0.0;
  bnrm2  = 0.0;
  error  = 0.0;
  rho    = 0.0;
  rho_1  = 0.0;

  dim    = length( b );
  x      = zeros(dim,1);
  p      = zeros(dim,1);
  q      = zeros(dim,1);
  r      = zeros(dim,1);
  z      = zeros(dim,1);

  % -----------------------------
  % Quick check of approximation.
  % -----------------------------

  bnrm2 = norm( b );
  if ( bnrm2 == 0.0 )
     bnrm2 = 1.0;
  end

  r = b;
  error = norm( r ) / bnrm2;
  if ( error < tol )
     return
  end

  % ----------------
  % Begin iteration.
  % ----------------

  if ( verbose )
     fprintf( '\n iter   prered    ||rres||    ||x||     norms \n\n' )
  end

  for iter = 1:max_it

     if ( verbose )
        if ( iter > 1 )
           fprintf( ' %4d  %+.2e  %.2e  %.2e  %.2e\n', iter, prered, error, norm(x), norms )
	else
           fprintf( ' %4d  %+.2e  %.2e  %.2e\n', iter, prered, error, norm(x) )
	end
     end

     %-----------------
     %  Preconditioning
     %-----------------

     if ( use_precond )
        if ( isa( precond, 'function_handle' ) )
           z = feval( precond, r );
        else
           z = precond * r;
	end
     else
        z = r;
     end
     rho = ( r' * z );

     % ---------------------
     % Reorthonormalization.
     % ---------------------

     if ( reorth )
        pnrz = sqrt( rho );
        if ( iter == 1 )
           R = [ z/pnrz ];
           Q = [ r/pnrz ];
	else
           z = z - R*(Q'*z);
           R = [ R, z/pnrz];
           Q = [ Q, r/pnrz];
        end
     end

     % -------------------------
     % Compute direction vector.
     % -------------------------

     if ( iter > 1 ),
        beta = rho / rho_1;
        p = z + beta * p;
     else
        p = z;
     end

     if ( isa( matvec, 'function_handle' ) )
        q = feval( matvec, p, S, Y );
     else
        q = matvec * p;       
     end
     pq    = p' * q;
     alpha = rho / pq;

     if ( TR_use || pq <= 0 )
        px  = p' * x;
        np2 = p' * p;
        nx2 = x' * x;
        alrad = (-px + sqrt( px^2 - np2 * ( nx2 - radius^2 ) ) ) / np2;
        if ( alpha > alrad || alpha <= 0 )
           at_boundary = 1;
           alpha       = alrad;
        else
           at_boundary = 0;
        end
     end
     
     % ---------------------
     % Update approximation.
     % ---------------------

     pr     = p' * r;
     x      = x + alpha * p;
     if ( verbose )
        norms = norm( alpha * p );
     end
     r      = r - alpha * q;
     prered = prered + alpha * pr - 0.5 * alpha^2 * pq;
     if ( prered < 0 || ~isreal( prered ) )
        disp( ' Problematic prered in pcg! ')
	prered
	alpha
	alrad
	if ( TR_use )
	   radius
	end
	pr
	pq
	normr = norm(r)
        keyboard%D
     end
     rnrm2  = norm( r );

     % ------------------
     % Check convergence.
     % ------------------

     error  = rnrm2 / bnrm2;
     if ( error <= tol || ( TR_use && at_boundary )   )
        break
     end
     rho_1 = rho;

  end

  % ------------------------
  % Final convergence check.
  % ------------------------

  if ( error > tol )
     if ( TR_use && at_boundary )
        if ( verbose )
           fprintf( ' %4d  %+.2e  %.2e  %.2e  %.2e  B\n\n', iter+1, prered, error, norm(x), norms )
        end
        flag = 2;
     else
        if ( verbose )
           fprintf( ' %4d  %+.2e  %.2e  %.2e  %.2e  M\n\n', iter+1, prered, error, norm(x), norms )
        end
        flag = 1;
     end
  else
     if ( verbose )
        fprintf( ' %4d  %+.2e  %.2e  %.2e  %.2e  C\n\n', iter+1, prered, error, norm(x), norms )
     end
  end

% --------
% End tcg.m 