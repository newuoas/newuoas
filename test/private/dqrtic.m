%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = dqrtic( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Variable dimension diagonal quadratic problem.
%
%   Source:  problem 157 (p.87) in
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
      if ( n < 1 )
         disp( [ ' ERROR in dqrtic: n = ', int2str(n),' should be > 0!' ] )
      end
   else
      n = 10;
   end
   varargout{1} = [ 2*ones( n, 1) ];            % x0
   varargout{2} = 0;                            % fstar
   varargout{3} = '';                           % xtype
   varargout{4} = [];                           % xlower
   varargout{5} = [];                           % xupper
   varargout{6} = [];                           % clower
   varargout{7} = [];                           % cupper
   varargout{8} = 'QUR2-AY-V-0';                % class

case 'objf'   % varargout = [ f, g, H ]

   x   = varargin{1};
   n   = length( x );
   for iel = 1:n
      eldom{ iel } = [ iel ];
   end
   varargout = opm_eval_cpsf( 'dqrtic', 'elobjf', x, eldom, nargout );

case 'elobjf' % varargout = [ fiel, giel, Hiel ]

   iel   = varargin{1};
   x     = varargin{2};
   riel  = x(1) - iel;
   varargout{1} = riel^2;
   if ( nargout > 1 )
      varargout{2} = 2 * riel;
%varargout{2}%D      
      if ( nargout > 2 )
         varargout{3} = 2;
      end
   end

return

end