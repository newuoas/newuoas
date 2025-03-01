%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = arglina( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Variable dimension full-rank linear least-squares problem.
%
%   Source: problem 32 in
%      J.J. More, B.S. Garbow and K.E. Hillstrom,
%      "Testing Unconstrained Optimization Software",
%      ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.
%   Also problem 80 in
%      A.R. Buckley,
%      "Test functions for unconstrained minimization",
%      TR 1989CS-3, Mathematics, statistics and computing centre,
%      Dalhousie University, Halifax (CDN), 1989.
%
%   If the dimension is unspecified, the default n = 10 is chosen.
%
%   Ph. Toint 17 VII 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   if ( length( varargin ) )
      n = varargin{1};
      if ( n < 1 )
         disp( [ ' ERROR in arglina: n = ', int2str(n),' but should satisfy n >= 1!' ] )
      end
   else
      n = 10;
   end
   m            = 2*n;
   varargout{1} = ones( n, 1 );                 % x0
   varargout{2} = m - n;                        % fstar
   varargout{3} = '';                           % xtype
   varargout{4} = [];                           % xlower
   varargout{5} = [];                           % xupper
   varargout{6} = [];                           % clower
   varargout{7} = [];                           % cupper
   varargout{8} = 'SUR2-AN-V-0';                % class

case 'objf'   % varargout = [ f, g, H ]

   x     = varargin{1};
   n     = length( x );
   m     = 2 * n;
   eldom = cell( m, 1 );
   for iel = 1:m
      eldom{ iel } = [ 1:n ];
   end
   varargout = opm_eval_cpsf( 'arglina', 'elobjf', x, eldom, nargout );

case 'elobjf' % varargout = [ fiel, giel, Hiel ]

   iel = varargin{1};
   x   = varargin{2};
   n   = length( x );
   m   = 2 * n;
   if ( iel <= n )
      riel = x(iel) - 1;
   else
      riel = -1;
   end
   for j = 1:n
      riel = riel - (2/m)*x(j);
   end
   varargout{1} = riel^2;
   if ( nargout > 1 )
      Jiel = zeros( n, 1 );
      if ( iel <= n )
         Jiel( iel ) = 1;
      end
      Jiel = Jiel - (2/m) * ones( n, 1 );
      varargout{2} = 2 * Jiel * riel;
      if ( nargout > 2 )
         varargout{3} = 2 * Jiel * Jiel.';
      end
   end

end

return

end