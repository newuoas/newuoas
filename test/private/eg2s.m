%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = eg2s( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   A variant on the EG2 problem.
%
%   Source:
%      ???
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
         disp( [ ' ERROR in eg2s: n = ', int2str(n),' should be > 0!' ] )
      end
   else
      n = 10;
   end
   varargout{1} = [ 8*ones( n, 1) ];            % x0
   varargout{2} = -8.29282661098;               % fstar
   varargout{3} = '';                           % xtype
   varargout{4} = [];                           % xlower
   varargout{5} = [];                           % xupper
   varargout{6} = [];                           % clower
   varargout{7} = [];                           % cupper
   varargout{8} = 'OUR2-AY-V-0';                % class

case 'objf'   % varargout = [ f, g, H ]

   x   = varargin{1};
   n   = length( x );
   for iel = 1:n-2
      eldom{ iel } = [ iel iel+1 iel+2 ];
   end
   varargout = opm_eval_cpsf( 'eg2s', 'elobjf', x, eldom, nargout, n );

case 'elobjf' % varargout = [ fiel, giel, Hiel ]

   iel = varargin{1};
   x   = varargin{2};
   n   = varargin{3};
   aux1  = sin( x(1) + x(2)^2 - 1 );
   varargout{1} = aux1 + (0.5/n) * sin(x(3)^2) ;
   if ( nargout > 1 )
      aux2 = cos( x(1) + x(2)^2 - 1 );
      varargout{2} = [ aux2; aux2*2*x(2); (0.5/n)*cos(x(3)^2)*2*x(3) ];
      if ( nargout > 2 )
         varargout{3} = [ -aux1       , -aux1*2*x(2)          , 0;
                          -aux1*2*x(2),  2*aux2-4*x(2)^2*aux1 , 0;
                              0,                  0           , (0.5/n)*(cos(x(3)^2)*2 -(2*x(3))^2*sin(x(3)^2))] ;
      end
   end

return

end