%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  varargout = cliff( action, varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Cliff problem in 2 variables.
%
%   Source: problem 206 in 
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

   varargout{1} = [ 0; -1];                     % x0
   varargout{2} = 0.19978661;                   % fstar
   varargout{3} = '';                           % xtype
   varargout{4} = [];                           % xlower
   varargout{5} = [];                           % xupper
   varargout{6} = [];                           % clower
   varargout{7} = [];                           % cupper
   varargout{8} = 'OUR2-AN-2-0';                % class

case 'objf'   % varargout = [ f, g, H ]

   x  = varargin{1};
   varargout = opm_eval_cpsf( 'cliff', 'elobjf', x, { [ 1 2 ] }, nargout );

case 'elobjf' % varargout = [ fiel, giel, Hiel ]

   iel  = varargin{1};
   x    = varargin{2};
   ee   = exp(20*(x(1)-x(2)));
   varargout{1} = ((x(1)-3)/100)^2 - (x(1) - x(2)) + ee;
   if ( nargout > 1 )
      varargout{2} = [ 0.0002*(x(1)-3)-1+20*ee; 1-20*ee ];
      if ( nargout > 2 )
	 varargout{3} = [  0.0002+400*ee, -400*ee;
	                   -400*ee       ,  400*ee ];
      end
   end
end

return

end