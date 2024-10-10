%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = eg2( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   A simple nonlinear problem given as an example in Section 1.2.4 of
%   the LANCELOT Manual.
%   The problem is non convex and has several local minima. The Hessian
%   is diagonal.
%
%   Source:
%      A.R. Conn, N. Gould and Ph.L. Toint,
%      "LANCELOT, A Fortran Package for Large-Scale Nonlinear Optimization
%      (Release A)"
%      Springer Verlag, 1992.
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
         disp( [ ' ERROR in eg2: n = ', int2str(n),' should be > 0!' ] )
      end
   else
      n = 10;
   end
   varargout{1} = [ 8*ones( n, 1) ];            % x0
   varargout{2} = -9.5;                         % fstar (n=10)
   varargout{3} = '';                           % xtype
   varargout{4} = [];                           % xlower
   varargout{5} = [];                           % xupper
   varargout{6} = [];                           % clower
   varargout{7} = [];                           % cupper
   varargout{8} = 'OUR2-AY-V-0';                % class

case 'objf'   % varargout = [ f, g, H ]

   x   = varargin{1};
   n   = length( x );
   for iel = 1:n
      eldom{ iel } = [ iel ];
   end
   varargout = opm_eval_cpsf( 'eg2', 'elobjf', x, eldom, nargout, n );

case 'elobjf' % varargout = [ fiel, giel, Hiel ]

   iel = varargin{1};
   x   = varargin{2};
   n   = varargin{3};
   if ( iel < n )
      aux1 = sin( x(1) + x(1)^2 - 1 );
      varargout{1} = aux1 ;
   else
      varargout{1} = 0.5 * sin( x(1)^2 );
   end
   if ( nargout > 1 )
      if ( iel < n )
         aux2 = cos( x(1) + x(1)^2 - 1 );
         varargout{2} = aux2 *( 1 +  2 * x(1) );
      else
         varargout{2} = cos(x(1)^2)*x(1);
      end
      if ( nargout > 2 )
         if ( iel < n )
            varargout{3} = -aux1-aux1*2*x(1) -aux1*2*x(1)+ 2*aux2-4*x(1)^2*aux1;
         else
            varargout{3} = .5*cos(x(1)^2)*2-.5*(2*x(1))^2*sin(x(1)^2);
         end
      end
   end

return

end