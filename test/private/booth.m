%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = booth( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Booth problem in 2 variables.
%   A simple quadratic in 2 variables.
%
%   Source: problem 36 in
%      A.R. Buckley,
%      "Test functions for unconstrained minimization",
%      TR 1989CS-3, Mathematics, statistics and computing centre,
%      Dalhousie University, Halifax (CDN), 1989.
%
%   Ph. Toint 17 VII 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   varargout{1} = [ 0; 0] ;                     % x0
   varargout{2} = 0;                            % fstar
   varargout{3} = '';                           % xtype
   varargout{4} = [];                           % xlower
   varargout{5} = [];                           % xupper
   varargout{6} = [];                           % clower
   varargout{7} = [];                           % cupper
   varargout{8} = 'QUR2-AN-2-0';                % class

case 'objf'   % varargout = [ f, g, H ]

   x  = varargin{1};
   varargout = opm_eval_cpsf( 'booth', 'elobjf', x, { [ 1 2 ] }, nargout );

case 'elobjf' % varargout = [ fiel, giel, Hiel ]

   iel  = varargin{1};
   x    = varargin{2};
   r1   =   x(1) + 2*x(2) - 7;
   r2   = 2*x(1) +   x(2) - 5;
   varargout{1} = r1^2 + r2^2;
   if ( nargout > 1 )
      J1 = [ 1; 2 ];
      J2 = [ 2; 1 ];
      varargout{2} = 2 * ( J1*r1 + J2*r2 );
      if ( nargout > 2 )
	    varargout{3} = 2 * ( J1*J1.' + J2*J2.' );
      end
   end
end

return

end