%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = bdarwhd( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   A variable-dimension quartic problem whose Hessian is an
%   arrow-head (downwards) with tridiagonal central part and border-width of 1.
%
%   Source: ???
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
      if ( n < 1 )
         disp( [ ' ERROR in bdarwhd: n = ', int2str(n),' but should satisfy n >= 1!' ] )
      end
   else
      n = 10;
   end
   varargout{1} = ones( n, 1 );                 % x0
   varargout{2} = 0;                            % fstar
   varargout{3} = '';                           % xtype
   varargout{4} = [];                           % xlower
   varargout{5} = [];                           % xupper
   varargout{6} = [];                           % clower
   varargout{7} = [];                           % cupper
   varargout{8} = 'OUR2-AY-V-0';                % class

case 'objf'   % varargout = [ f, g, H ]

   x     = varargin{1};
   n     = length( x );
   eldom = cell( n-2, 1 );
   for iel = 1:n-2
      eldom{ iel } = [ iel iel+1 n ];
   end
   varargout = opm_eval_cpsf( 'bdarwhd', 'elobjf', x, eldom, nargout );

case 'elobjf' % varargout = [ fiel, giel, Hiel ]

   iel  = varargin{1};
   x    = varargin{2};
   riel = x(1) + x(2) + x(3);
   varargout{1} = riel^4;
   if ( nargout > 1 )
      varargout{2} = 4 * riel^3 * ones( 3, 1 ); 
      if ( nargout > 2 )
	 varargout{3} = 12 * riel^2 * ones( 3, 3 );
      end
   end
   
end

return

end