function [ xs, barfs, bargs, it, omegas ] = ilmqn_opm( prob, x0, Delta, epsilon, maxit, inexact )

%  Solves the unconstrained problem min_x f(x) starting from x0 using a variable accuracy strategy.

%  Ph. Toint 3 V 2018

verbose   = 0;

n         = length( x0 );    % problem size
pmax      = min( 5, n );     % max number of secant pairs-1
%update    = 'BFGS';
%update   = 'SR1';
update    = 'SR1TR';
globalize = 'TR';
%globalize = 'LS';

S       = [];              % S and Y gather the secant pairs for preconditioning
Y       = [];
tau     = [];
it      = 0;               % outer iteration counter
x       = x0;
Delprev = Delta;
Delmax  = 1e8;
eta1    = 0.25;
eta2    = 0.95;
eta0    = 0.4 * eta1;
kappa_g = 0.5 * ( 1 - eta2 );
gamma1  = 1/3;
gamma2  = 3;
itmdmax = 10; %% ???
itouter = 0;
xs      = [ x0 ];

%Delta = 1;   %%%%  IMPROVE !!!

%  Compute approximate f and g up to accuracy omega_f and omega_g, respectively.

if ( inexact )
   omega_f  = 0.01;
   omega_g  = kappa_g;
else
   omega_f  = 0;
   omega_g  = 0;
end
[ barf, barg, eff_omega_f ] = iobjf_opm( prob, x, omega_f, omega_g );
nf          = 1;
ng          = 1;
nbarg       = norm( barg );
relf        = max( 1, abs( barf ) );
barfs       = [ barf ];
bargs       = [ barg ];
omegas(1,1) = omega_f;
omegas(2,1) = omega_g;
preredp     = 1;

[ fx, gx ] = feval( prob, 'objf', x );
aerrf  = abs( fx - barf );
rerrg  = norm( gx - barg ) / max( 1, nbarg );
barfs  = [ barfs barf ];
bargs  = [ bargs barg ];

M = eye(n);%D

%  Initial printout

if ( verbose )
   switch ( globalize )
   case 'TR'
      fprintf( '\n it   itmd    nf    ng      barf       aerrf     ||barg||    rerrg   omega_f    prered     ||s||      Delta        rho \n\n' )
   case 'LS'
      fprintf( '\n it   itmd    nf    ng      barf       aerrf     ||barg||    rerrg   omega_f    prered       ||s||     stepsize \n\n' )
   end
end

%  Loop for outer iterations

while( it <= maxit )

%  Printout

   if ( verbose )
      if ( it )
         switch ( globalize )
	 case 'TR'
            fprintf( '%4d  %4d  %4d  %4d  %+.5e  %.1e  %.5e  %.1e  %.1e  %.3e  %.3e  %.3e  %+.3e\n',  ...
                   it, itmd, nf, ng, barf, aerrf, nbarg, rerrg, omega_f, prered, norms, Delprev, rho );
         case 'LS'
            fprintf( '%4d  %4d  %4d  %4d  %+.5e  %.1e  %.5e  %.1e  %.1e  %.3e  %.3e  %.3e\n',         ...
                   it, itmd, nf, ng, barf, aerrf, nbarg, rerrg, omega_f, prered, norms, stepsize );
         end	  
      else
         fprintf( '%4d        %4d  %4d  %+.5e  %.1e  %.5e  %.1e  %.1e\n',                             ...
                it, nf, ng, barf, aerrf, nbarg, rerrg, omega_f);
      end
   end

   %  Possibly terminate

   if ( ( 1 + omega_g ) * nbarg <= epsilon )
      if ( verbose )
         fprintf( '\n  Done.\n' )
      end
      return;
   end
   it = it + 1;

   %  Compute the step and the model decrease.

   itmd   = 0;
   while ( itmd <= itmdmax )

      itmd = itmd + 1;

%     Compute a search direction by applying a QN(Y,S) method

      switch ( update )
      case 'BFGS'
         s     = lmtrs( barg, Y, S, Delta, 50, 0.001 );
         norms = norm( s );
         switch ( globalize )
	 case 'LS'
            barfplus = iobjf_opm( prob, x+s, omega_f );
	    nf       = nf + 1;
	    sbarg    = s' * barg;
	    stepsize = -0.5 * sbarg / ( barfplus - barf - sbarg );
            s        = stepsize * s;
            norms    = stepsize * norms;
	 case 'TR'
         end
         prered = -0.5 * barg' * s;
      case 'SR1'
         if ( size( S, 2 ) )
            gamma = 2*(S(:,1)'*Y(:,1))/(S(:,1)'*S(:,1));
            gamma = 0.01;
            M     = S' * Y;
            Psi   = Y - gamma * S;
	    if ( rank( Psi ) ~=  size(S,2) )
	       rank(Psi)
	       rank(S)
	       size(S,2)
	       pause
	    end
            invM  = tril(M) + tril(M,-1)' - gamma * S' * S;
            [ lambda, s ] = obs( barg, [], [], Delta, gamma, Psi, invM );
            norms = norm( s );
            prered = -0.5 * barg' * s;
	 else
	    s      = -barg * Delta / nbarg;
            norms  = Delta;
            prered = -barg' * s;
	 end
      case 'SR1TR'
         if ( size( S, 2 ) )
	 
NEWTON = 1;
if ( NEWTON )
    [ ff, gg, H ] = feval( prob, 'objf', x );
    [ s, ~, ~, ~, prered ] = pcg( H, -barg, 0, 3*n, 0.0001, Delta, 0, S, Y );
else
            [ s, ~, ~, ~, prered ] = pcg( @prodSR1, -barg, 0, 3*n, 0.0001, Delta, 0, S, Y );
end%
            norms    = norm( s );
	    sbarg    = barg'*s;
	    stepsize = 0.5 * sbarg / ( prered + sbarg );
	    if ( stepsize > 0 && stepsize < 1 )
	       s     = stepsize * s;
	       norms = stepsize * norms;
	    end
	 else
	    s      = -barg * Delta / nbarg;
            norms  =  Delta;
            prered = -barg' * s;
	 end
      end
      
%     Exit of the model decrease loop if prered is large enough

      mu = omega_f / prered;
      if ( mu <= eta0 )
         break
      end

      %  Reduce omega.
      
      omega_f = min( [ nbarg^2, 0.5 * omega_f, 0.5*eta0*prered ] );

      %  Recompute the current barf and barg.

      [ barf, barg, eff_omega_f ] = iobjf_opm( prob, x, omega_f, omega_g );
      nbarg = norm( barg );
      nf    = nf + 1;
      ng    = ng + 1;
      
   end

   %  Define the next omega and compute associated f and g.

   %  omega must be small enough compared to prered and, at
   %  the same time, small enough to ensure the suitable
   %  relative errors on f and g.  The present version is
   %  somewhat heuristic. The geometric mean between prered and preredp
   %  ensures that omega does not grow too fast, leading to oscillations

   xnew = x + s;
   if ( inexact )
      if ( itmd < 2 )
          omeganew = min( [ 0.1,  eta0*sqrt(preredp*prered), nbarg^2 ] );
      else
          omeganew = min( [ omega_f, eta0*sqrt(preredp*prered), nbarg^2 ] );
      end
   else
      omeganew = 0;
   end
   [ barfnew, bargnew, eff_omega_f ] = iobjf_opm( prob, xnew, omeganew, omega_g ); 
   nf = nf + 1;
   ng = ng + 1;
   y = bargnew - barg;
   nbargnew = norm( bargnew );
   [ fxnew, gxnew ] = feval( prob, 'objf', xnew );
   aerrfnew = abs( fxnew - barfnew );
   rerrgnew = norm( gxnew - bargnew ) / max( 1, nbargnew);
   barfs    = [ barfs barfnew ];
   bargs    = [ bargs bargnew ];
   xs       = [ xs xnew ];
   
   %  Update the trust-region radius and accept/reject the trial point.

   switch( globalize )
   case 'LS'
      x       = xnew;
      barf    = barfnew;
      barg    = bargnew;
      nbarg   = nbargnew;
      aerrf   = aerrfnew;
      rerrg   = rerrgnew;
      fx      = fxnew;
      gx      = gxnew;
      omega_f = omeganew;
      rho     = 1;         % a convention
   case 'TR'
      ared    = barf - barfnew;
      rho     = (ared + 1.e-16) / (prered + 1.e-16);
      Delprev = Delta;

      if ( rho > eta1 )             %  Successful iteration

         x       = xnew;
         barf    = barfnew;
         barg    = bargnew;
         nbarg   = nbargnew;
         aerrf   = aerrfnew;
         rerrg   = rerrgnew;
         fx      = fxnew;
         gx      = gxnew;
         omega_f = omeganew;

         if ( rho > eta2 )          %  Very successful iteration
            Delta = min( gamma2 * norms, Delmax );
         end
	 
      else                          %  Unsuccessful iteration
      
         Delta   = gamma1 * min( norms, Delta );
         omega_f = 0.5 * omeganew;  % a precaution, not strictly necessary
      end
      
   end
   
   %  Accumulate the secant pairs

%   y    = bargnew - barg;
   skip = ( norm( y ) > 1e10 * max(1,norm( Y ))) || ( strcmp( update, 'BFGS' ) && s' * y  < 1e-10 * norms^2 );
   if ( ~skip )
      S(1:n,end+1) = s;
      Y(1:n,end+1) = y;
      tau(end+1)   = omeganew;
      if ( size( S, 2 ) > pmax )
         S   = S(1:n,end-pmax:end );
         Y   = Y(1:n,end-pmax:end );
         tau = tau(end-pmax:end );
      end
   end

   %  Update omegas (a bit heuristic for now).
   
   omega_g = min( [ kappa_g, 2*omega_f/(nbarg*norms), omega_f^(2/3) ] );

   preredp = prered;
   
   omegas(1,end+1) = omega_f;
   omegas(2,end+1) = omega_g;

end

return

end
