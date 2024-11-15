%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = zangwil2( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   The Zangwill quadratic problem is 2 variables.
%
%   Source: problem 7 (p. 102) in 
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

   varargout{1} = [ 3; 8];                      % x0
   varargout{2} = -18.2;                        % fstar
   varargout{3} = '';                           % xtype
   varargout{4} = [];                           % xlower
   varargout{5} = [];                           % xupper
   varargout{6} = [];                           % clower
   varargout{7} = [];                           % cupper
   varargout{8} = 'QUR2-AN-2-0';                % class

case 'objf'   % varargout = [ f, g, H ]

   x  = varargin{1};
   varargout = opm_eval_cpsf( 'zangwil2', 'elobjf', x, { [ 1 2 ] }, nargout );

case 'elobjf' % varargout = [ fiel, giel, Hiel ]

   x    = varargin{2};
   varargout{1} = (1/15)*( 16*x(1)^2 + 16*x(2)^2 -8*x(1)*x(2) -56*x(1) -256*x(2) + 991 );
   if ( nargout > 1 )
      varargout{2} = (1/15) * [ 32*x(1)-8*x(2)-56; 32*x(2)-8*x(1)-256 ];
      if ( nargout > 2 )
	 varargout{3} = (1/15) * [ 32, -8;
	                           -8, 32 ];
      end
   end
end

return

end