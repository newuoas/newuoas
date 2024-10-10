function [ err_gf, err_Hf, err_gc, err_Hc, objflin, conslin ] = opm_test ( probname, verbose, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Test the derivatives of the functions for problem probname, using finite differences in continuous
%   varaiables at a random point feasible for the problems bounds. Also determines if the objective
%   function and constraints are (likely to be) linear or not by verifying if the Hessian is zero
%   at the (random) evaluation points.
%
%   INPUT:
%
%   probname: the name of the considered optimization problem
%   verbose : 1 if output and the standard output is desired, 0 for silent execution
%   varargin: possible parameters.  If specified varargin{1} must be the problem's dimension.
%
%   OUTPUT:
%
%   err_gf : the max absolute value of the error on the objective function's gradient components
%   err_Hf : the max absolute value of the error on the objective function's Hessian components
%   err_gc : the max absolute value of the error on the constraints' function's gradient components
%   err_Hc : the max absolute value of the error on the constraints' function's Hessian components
%   objflin: 1 if the objective function appears to be linear, 0 otherwise
%   conslin: 1 if all constraint functions appears to be linear, 0 otherwise
%
%   Programming : Ph. Toint and S. Gratton, July 2018
%   This version: 17 VII 2018
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ep = 1e-12;

%  Get problem data.

if ( nargin > 2 )
   [ x0, fstar, xtype, xlower, xupper, clower, cupper, class ] = feval( probname, 'setup', varargin{:} );
else
   [ x0, fstar, xtype, xlower, xupper, clower, cupper, class ] = feval( probname, 'setup' );
end
n  = length( x0 );
if ( length( xtype ) )
   xcont = find( xtype == 'c' );
else
   xcont = [1:n];
end
%xnotcont = setdiff( [1:n], xcont );
ncont    = length( xcont );

%  Choose a random bound-feasible evaluation point.

x  = rand( n, 1 );

if ( length( xlower ) )
   x  = min( [ xupper';  max( [ x'; xlower' ] ) ] )';
end

%%%%%%%%%%%%%%%%%%%%%%%%% Objective function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   Get the values of the objective function at x.

[ f, g, H ] = feval( probname, 'objf', x );

%   Test its gradient.

v       = zeros( n, 1 );
err_gf  = 0;
err_gfi = zeros( n, 1 );
for i   = 1:n
   if ( length( xlower ) == 0 || abs( xupper( i ) - xlower( i ) ) > 1e-15 )
      v( i ) = 1;
      [ fp, gp, Hp ] = feval( probname, 'objf', x + ep*sqrt(-1)*v );
      err_gfi( i )   = abs( (1/ep)*imag( fp ) - g(i) );
      v( i ) = 0;
   end
end

if ( verbose )
   disp(   ' ' )
   disp(  [ '*** Testing problem ', probname, ' ***' ] )
   disp(    ' ' )
   disp(    '    objective function : errors on the gradient' )
   disp(    ' ' )
   fprintf( '                 ' ); fprintf( '%0.3e  ', err_gfi(1:n) )
   disp( ' ')
end
err_gf = norm( err_gfi, Inf );

%   Test its Hessian.

err_Hf = zeros( n, n );
for i = 1:n
   if ( length( xlower ) == 0 || abs( xupper( i ) - xlower( i ) ) > 1e-15 )
      for j = 1:n
         if ( length( xlower ) == 0 || abs( xupper( j ) - xlower( j ) ) > 1e-15 )
            v( j ) = 1;
            [ fp, gp, Hp ] = feval( probname, 'objf', x + ep*sqrt(-1)*v );
            err_Hf( i, j ) = abs( (1/ep)*imag( gp(i) ) - H(i,j) );
            v( j ) = 0;
	 end
      end
   end
end
if ( verbose )
   disp(    ' ' )
   disp(    '    objective function : errors on the Hessian' )
   disp(    ' ' )
   for i = 1:n
   fprintf( '      row %2d  :  ', i ); fprintf( '%0.3e  ', err_Hf(i,1:n) ); fprintf( '\n' )
   end
end
err_Hf = norm( err_Hf, Inf );

%   Hopefully detect linearity.

objflin =( norm( H, Inf ) <= n * ep && norm( Hp, Inf ) <= n * ep );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Constraints  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ncon    = length( clower );
conslin = 1;
err_gc  = 0;
err_Hc  = 0;

for icon = 1:ncon

%  Get the values of the icon-th constraint at x.

   [ c, gc, Hc ]  = feval( probname, 'consf', icon, x );

%  Test its gradient.

   err_gci = zeros( n, 1 );
   for i  = 1:n
      xlower
      if ( length( xlower ) == 0 || abs( xupper( i ) - xlower( i ) ) > 1e-15 )
         v( i ) = 1;
         [ fp, gp, Hp ] = feval( probname, 'consf', icon, x + ep*sqrt(-1)*v );
         err_gci( i )   = abs( (1/ep)*imag( fp ) - gc(i) );
         v( i ) = 0;
      end
   end
   err_gc = max( err_gc, err_gci );

   if ( verbose )
      disp(    ' ' )
      disp(    [ '    constraint , int2str( icon ),': errors on the gradient' ])
      disp(    ' ' )
      fprintf( '                 ' ); fprintf( '%0.3e  ', err_gci(1:n) )
      disp(    ' ' )
   end

%  Test its Hessian.

   err_Hci = zeros( n, n );
   for i = 1:n
      if ( length( xlower ) == 0 || abs( xupper( i ) - xlower( i ) ) > 1e-15 )
         for j = 1:n
            if ( length( xlower ) == 0 || abs( xupper( j ) - xlower( j ) ) > 1e-15 )
               v( j ) = 1;
               [ fp, gp, Hp ]  = feval( probname, 'consf', icon, x + ep*sqrt(-1)*v );
               err_Hci( i, j ) = abs( (1/ep)*imag( gp(i) ) - Hc(i,j) );
               v( j ) = 0;
	    end
         end
      end
   end
   if ( verbose )
      disp(    ' ' )
      disp(    '    constraint', int2str(icon),': errors on the Hessian' )
      disp(    ' ' )
      for i = 1:n
         fprintf( '      row %2d  :  ', i ); fprintf( '%0.3e  ', err_Hci(i,1:n) ); fprintf( '\n' )
      end
   end
   err_Hc = max( err_Hc, norm( err_Hci, Inf ) );

   %   Hopefully detect linearity.

   conslin = conslin && ( norm( H, Inf ) <= n * ep && norm( Hp, Inf ) <= n * ep );

end

if ( verbose )
   disp( ' ------------------' )
end
