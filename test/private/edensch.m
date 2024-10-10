%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = edensch( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   The extended Dennis-Schnabel problem, as defined by Li.
%
%   Source:
%      "The secant/finite difference algorithm for solving sparse
%      nonlinear systems of equations",
%      SIAM Journal on Optimization, (to appear), 1990.
%
%   If the dimension is unspecified, the default n = 10 is chosen.
%
%   Ph. Toint and S. Gratton, 22 VII 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   if ( length( varargin ) )
      n = varargin{1};
      if ( n < 1 )
         disp( [ ' ERROR in edensch: n = ', int2str(n),' should be > 0!' ] )
      end
   else
      n = 5;
   end
   varargout{1} = [ 8*ones( n, 1) ];            % x0
   varargout{2} = 33.296001998982;              % fstar (n=5)
   varargout{3} = '';                           % xtype
   varargout{4} = [];                           % xlower
   varargout{5} = [];                           % xupper
   varargout{6} = [];                           % clower
   varargout{7} = [];                           % cupper
   varargout{8} = 'OUR2-AY-V-0';                % class

case 'objf'   % varargout = [ f, g, H ]

   x   = varargin{1};
   n   = length( x );
   for iel = 1:n-1
      eldom{ iel } = [ iel iel+1 ];
   end
   eldom{n} = n;
   varargout = opm_eval_cpsf( 'edensch', 'elobjf', x, eldom, nargout, n );

case 'elobjf' % varargout = [ fiel, giel, Hiel ]

   iel = varargin{1};
   x   = varargin{2};
   n   = varargin{3};
   if ( iel < n )
      aux = x(1) * x(2) - 2 * x(2);
      varargout{1}  = ( x(1)- 2 )^4 + aux^2 + ( x(2) + 1 )^2;
      if ( nargout > 1 )
         varargout{2} = [4*(x(1)- 2)^3+2*aux*x(2); 2*aux*(x(1)-2)+2*(x(2)+1)];
         if ( nargout > 2 )
            varargout{3} = [12*(x(1)- 2)^2+2*x(2)^2, 2*aux+2*(x(1)-2)*x(2);
                             2*aux+2*(x(1)-2)*x(2)  , 2*(x(1)-2)*(x(1)-2)+2 ];
         end
      end
   else
      varargout{1} = (-2)^4;        % (x(1)-2)^4;
      if ( nargout > 1 )
         varargout{2} = 0;          % 4*(x(1)-1)^3;
	 if ( nargout > 2 )
	    varargout{3} = 0;       % 12*(x(1)-1)^2;
	 end
      end
   end
return

end