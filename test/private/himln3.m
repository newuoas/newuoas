%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = himln3( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   A small 2D problem.
%
%   Source: problem 8 (p. 60) in
%      A.R. Buckley,
%      "Test functions for unconstrained minimization",
%      TR 1989CS-3, Mathematics, statistics and computing centre,
%      Dalhousie University, Halifax (CDN), 1989.
%
%   Ph. Toint 20 VII 2018.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch ( action )

case 'setup' % varargout = [ x0, fstar, xtype, xlower, xupper, clower, cupper, class, error ]

   varargout{1} = [ 0; 2 ];                     % x0
   varargout{2} = -1;                           % fstar
   varargout{4} = [];                           % xlower
   varargout{5} = [];                           % xupper
   varargout{6} = [];                           % clower
   varargout{7} = [];                           % cupper
   varargout{8} = 'OUR2-AN-2-0';                % class

case 'objf'   % varargout = [ f, g, H ]

   x   = varargin{1};
   varargout = opm_eval_cpsf( 'himln3', 'elobjf', x, { [ 1 2 ] }, nargout);

case 'elobjf' % varargout = [ fiel, giel, Hiel ]

   x     = varargin{2};
   varargout{1} = x(1)^3 + x(2)^2 - 3*x(1) - 2*x(2) + 2;
   if ( nargout > 1 )
      varargout{2} = [ 3*x(1)^2-3; 2*x(2)-2 ];
      if ( nargout > 2 )
         varargout{3} = [ 6*x(1), 0;
	                   0    , 2 ];
      end
   end
end

return

end