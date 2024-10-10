%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = broyden3d( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Broyden tridiagonal problem in variable dimension.  This is a nonlinear
%   least-squares problem with n groups.
%
%   Source: Problem 33 in
%      J.J. More', B.S. Garbow and K.E. Hillstrom,
%      "Testing Unconstrained Optimization Software",
%      ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.
%   Also problem 78 in 
%      A.R. Buckley,
%      "Test functions for unconstrained minimization",
%      TR 1989CS-3, Mathematics, statistics and computing centre,
%      Dalhousie University, Halifax (CDN), 1989.
%
%   If the dimension is unspecified, the default n = 10 is chosen.
%
%   Ph. Toint 22 VII 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   if ( length( varargin ) )
      n = varargin{1};
      if ( n < 3 )
         disp( [ ' ERROR in broyden3d: n = ', int2str(n),' but should satisfy n >= 3!' ] )
      end
   else
      n = 10;
   end
   varargout{1} = [ 0; ones( n-2, 1 ); 0 ];     % x0
   varargout{2} = 0.93929749097;                % fstar (n=10) ????
   varargout{3} = '';                           % xtype
   varargout{4} = [ 0; -Inf*ones( n-2, 1 ); 0 ];% xlower
   varargout{5} = [ 0;  Inf*ones( n-2, 1 ); 0 ];% xupper
   varargout{6} = [];                           % clower
   varargout{7} = [];                           % cupper
   varargout{8} = 'SXR2-AY-V-0';                % class

case 'objf'   % varargout = [ f, g, H ]

   x       = varargin{1};
   n       = length( x );
   eldom   = cell( n-2, 1 );
   for iel = 1:n-2
      eldom{ iel } = [ iel iel+1 iel+2 ];
   end
   varargout = opm_eval_cpsf( 'broyden3d', 'elobjf', x, eldom, nargout );

   %  Zero the derivatives of fixed variables.
   
   if ( nargout > 1 )
      varargout{2}(1) = 0;
      varargout{2}(n) = 0;
      if ( nargout > 2 )
         varargout{3}( [ 1, n ], : ) = zeros( 2, n );
         varargout{3}( :, [ 1, n ] ) = zeros( n, 2 );
      end
   end
   
case 'elobjf' % varargout = [ fiel, giel, Hiel ]

%  iel     = varargin{1};
   x       = varargin{2};
   riel    = ( 3 - 2*x(2) )*x(2) - x(1) - 2*x(3) + 1;
   varargout{1} = riel^2;
   if ( nargout > 1 )
      Jiel  = [ -1; 3-4*x(2); -2 ];
      varargout{2} = 2 * Jiel * riel;
      if ( nargout > 2 )
         Hiel = [ 0,  0, 0;
                  0, -4, 0;
	  	  0,  0, 0 ];
         varargout{3} = 2 * ( Jiel*Jiel.' + riel*Hiel );
      end
   end
end

return

end