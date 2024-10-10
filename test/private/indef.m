%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = indef( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   A nonconvex problem which has an indefinite Hessian at
%   the starting point.
%
%   Source:
%      Nick Gould, CUTE, Oct 1992.
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
      if ( n < 2 )
         disp( [ ' ERROR in indef: n = ', int2str(n),' does not satisfy n > 1!' ] )
      end
   else
      n = 5;
   end
   varargout{1} = [ 1:n ]'/(n+1);               % x0
   varargout{2} = -500.4416811112;             % fstar (n=5)
   varargout{3} = '';                           % xtype
   varargout{4} = [];                           % xlower
   varargout{5} = [];                           % xupper
   varargout{6} = [];                           % clower
   varargout{7} = [];                           % cupper
   varargout{8} = 'OUR2-AY-V-0';                % class

case 'objf'   % varargout = [ f, g, H ]

   x  = varargin{1};
   n  = length( x );
   for iel = 1:n                         %
      if ( iel == 1 || iel == n )
         eldom{ iel } = [ iel ];
      else
         eldom{ iel } = [ 1 iel n ];
      end
   end
   varargout = opm_eval_cpsf( 'indef', 'elobjf', x, eldom, nargout, n );

case 'elobjf' % varargout = [ fiel, giel, Hiel ]

   iel = varargin{1};
   x   = varargin{2};
   n   = varargin{3};
   lx  = length( x );
   scale = 100;
   if ( iel == 1 || iel == n )
      varargout{1} = scale * sin( x(1)/scale ) ;
   else
      aux1 = cos( -x(1)+2*x(2)-x(3) );
      varargout{1} = 0.5 * aux1 + scale * sin( x(2)/scale ) ;
   end
   if ( nargout > 1 )
      if ( iel == 1 || iel == n )
         varargout{2} = cos( x(1)/scale ) ;
      else
         aux2 = sin( -x(1)+2*x(2)-x(3) );
         varargout{2} = [ 0.5*aux2; -aux2 + cos( x(2)/scale ); 0.5*aux2 ];
      end
      if ( nargout > 2 )
         if ( iel == 1 || iel == n )
            varargout{3} = -(1/scale)*sin( x(1)/scale );
         else
            varargout{3} = [ -0.5*aux1,    aux1                             , -0.5*aux1;
                                  aux1, -2*aux1-(1/scale)*sin( x(2)/scale ) ,      aux1;
                             -0.5*aux1,    aux1                             , -0.5*aux1 ];
         end
      end
   end

end

return

end