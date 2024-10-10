%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = rosenbr( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   The ever famous 2 variables Rosenbrock "banana valley" problem
%
%   Source: problem 1 in
%      J.J. More, B.S. Garbow and K.E. Hillstrom,
%      "Testing Unconstrained Optimization Software",
%      ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.
%   Also problem 1 (p. 89) in
%      A.R. Buckley,
%      "Test functions for unconstrained minimization",
%      TR 1989CS-3, Mathematics, statistics and computing centre,
%      Dalhousie University, Halifax (CDN), 1989.
%
%   If the dimension is unspecified, the default n = 10 is chosen.
%
%   Ph. Toint 26 VII 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   if ( length( varargin ) )
      n = varargin{1};
      if ( n < 2 )
         disp( [ ' ERROR in rosenbr: n = ', int2str(n), ' must be > 1!' ] )
      end
   else
      n = 10;
   end
   varargout{1} = [ -1.2; ones( n-1,1) ];       % x0
   varargout{2} = 0;                            % fstar
   varargout{3} = '';                           % xtype
   varargout{4} = [];                           % xlower
   varargout{5} = [];                           % xupper
   varargout{6} = [];                           % clower
   varargout{7} = [];                           % cupper
   varargout{8} = 'SUR2-AY-V-0';                % class

case 'objf'   % varargout = [ f, g, H ]

   x     = varargin{1};
   n     = length( x );
   for iel = 1:n-1
      eldom{ iel } = [ iel iel+1 ];
   end
   varargout = opm_eval_cpsf( 'rosenbr', 'elobjf', x, eldom, nargout, n );

case 'elobjf' % varargout = [ fiel, giel, Hiel ]

   iel  = varargin{1};
   x    = varargin{2};
   r1   = 10*(x(2)-x(1)^2);
   r2   = 1-x(1);
   varargout{1} = r1^2 + r2^2;
   if ( nargout > 1 )
      J1 = 10*[-2*x(1); 1 ];
      J2 = [ -1; 0 ];
      varargout{2} = 2* ( J1*r1 + J2*r2 );
      if ( nargout > 2 )
         H1 = [ -20, 0; 0, 0  ];
	 varargout{3} = 2 * ( J1*J1.' + r1*H1 + J2*J2.' );
      end
   end
end

return

end