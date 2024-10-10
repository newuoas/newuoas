%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = freuroth( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   The extended Freudentstein and Roth test problem.
%
%   Source: problem 2 in
%      J.J. More', B.S. Garbow and K.E. Hillstrom,
%      "Testing Unconstrained Optimization Software",
%      ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.
%   Also problem 24 in
%      A.R. Buckley,
%      "Test functions for unconstrained minimization",
%      TR 1989CS-3, Mathematics, statistics and computing centre,
%      Dalhousie University, Halifax (CDN), 1989.
%   and problem 33 in
%      Ph.L. Toint,
%      "Test problems for partially separable optimization and results
%      for the routine PSPMIN",
%      Report 83/4, Department of Mathematics, FUNDP (Namur, B), 1983.
%
%   If the dimension is unspecified, the default n = 10 is chosen.
%
%   Ph. Toint and S. Gratton,  22 VII 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   if ( length( varargin ) )
      n = varargin{1};
      if ( n < 2 )
         disp( [ ' ERROR in freuroth: n = ', int2str(n),' should be > 1!' ] )
      end
   else
      n = 4;
   end
   varargout{1} = -2*ones( n, 1);               % x0
   varargout{2} = 284.044708637;                % fstar (n=4)
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
   varargout = opm_eval_cpsf( 'freuroth', 'elobjf', x, eldom, nargout, n );

case 'elobjf' % varargout = [ fiel, giel, Hiel ]

   iel = varargin{1};
   x   = varargin{2};
   n   = varargin{3};
   r1 = x(1) - 13 + ((5-x(2))*x(2)-2)*x(2);
   r2 = x(1) - 29 + ((x(2)+1)*x(2)-14)*x(2);
   varargout{1} = r1^2 + r2^2;
   if ( nargout > 1 )
      J1 = [ 1; 10*x(2)-3*x(2)^2-2 ];
      J2 = [ 1; 3*x(2)^2+2*x(2)-14 ];
      varargout{2} = 2 *( J1 * r1 + J2 * r2 );
      if ( nargout > 2 )
         H1 = [ 0,        0 ;
                0, 10-6*x(2) ];
         H2 = [ 0,        0 ;
                0,  6*x(2)+2 ];
         varargout{3} = 2 * ( J1*J1.' + r1*H1 + J2*J2.' + r2*H2 );
      end
   end

return

end